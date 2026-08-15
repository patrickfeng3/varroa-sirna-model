#!/usr/bin/env python3
"""Canonical Stage 03 official-stepRNA duplex-geometry reconstruction."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Iterable, Iterator


STRANDS = ("sense", "antisense")
EXPECTED_CATEGORIES = {
    "mapping_mode": {"exact", "1mm"},
    "virus_assignment": {"assigned", "ambiguous_multi_virus"},
    "strand": {"sense", "antisense", "ambiguous"},
}
REQUIRED_OFFICIAL_SUFFIXES = (
    "_overhang.csv",
    "_unique_overhang.csv",
    "_overhang_type.csv",
    "_passenger_length.csv",
    "_passenger_number.csv",
)


class Stage03Error(RuntimeError):
    """Structured Stage 03 execution or parsing failure."""


@dataclass(frozen=True)
class Stage03Config:
    focal_lengths: tuple[int, ...]
    passenger_length_min: int
    passenger_length_max: int
    joint_5p_distance: int
    joint_3p_distance: int
    required_steprna_version: str


def load_config(path: Path) -> Stage03Config:
    data = json.loads(path.read_text())["stage03"]
    config = Stage03Config(
        focal_lengths=tuple(int(x) for x in data["focal_lengths"]),
        passenger_length_min=int(data["passenger_length_min"]),
        passenger_length_max=int(data["passenger_length_max"]),
        joint_5p_distance=int(data["joint_5p_distance"]),
        joint_3p_distance=int(data["joint_3p_distance"]),
        required_steprna_version=str(data["steprna_version"]),
    )
    if config.focal_lengths != (23, 24):
        raise ValueError("Stage 03 canonical focal lengths must be [23, 24]")
    if (config.passenger_length_min, config.passenger_length_max) != (15, 30):
        raise ValueError("Stage 03 canonical passenger range must be 15-30")
    if (config.joint_5p_distance, config.joint_3p_distance) != (2, -2):
        raise ValueError("Stage 03 canonical joint geometry must be (+2,-2)")
    if config.required_steprna_version != "1.0.6":
        raise ValueError("Stage 03 requires official stepRNA 1.0.6")
    return config


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def opposite_strand(strand: str) -> str:
    if strand == "sense":
        return "antisense"
    if strand == "antisense":
        return "sense"
    raise ValueError(f"invalid strand: {strand}")


def stable_identifier(prefix: str, *parts: object) -> str:
    digest = hashlib.sha256("\x1f".join(str(x) for x in parts).encode()).hexdigest()[:24]
    return f"{prefix}_{digest}"


def run_identifier(sample: str, analysis_unit: str, length: int, strand: str) -> str:
    clean = lambda value: re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    suffix = "S" if strand == "sense" else "AS"
    return f"{clean(sample)}__{clean(analysis_unit)}__{length}{suffix}"


def safe_fraction(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        if value.is_integer():
            return str(int(value))
        return format(value, ".12g")
    return str(value)


def write_table(path: Path, rows: list[dict[str, object]], fields: list[str], gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    opener = gzip.open if gzip_output else open
    kwargs = {"mode": "wt", "newline": ""} if gzip_output else {"mode": "w", "newline": ""}
    with opener(temporary, **kwargs) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")
    os.replace(temporary, path)


def read_eligibility(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_feature_rows(root: Path, samples: Iterable[str]) -> Iterator[dict[str, str]]:
    for sample in sorted(samples):
        path = root / "tables" / sample / f"{sample}.read_level_features.tsv.gz"
        with gzip.open(path, "rt", newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def collapse_inputs(
    eligibility: list[dict[str, str]],
    feature_rows: Iterable[dict[str, str]],
    config: Stage03Config,
) -> dict[str, object]:
    metadata = {
        (row["sample"], row["analysis_unit"]): row
        for row in eligibility if is_true(row.get("primary_eligible"))
    }
    primary_pairs = set(metadata)
    sequence_abundance: dict[tuple[str, str, str, int], Counter[str]] = defaultdict(Counter)
    categories: dict[str, set[str]] = defaultdict(set)
    unexpected_sequences: set[str] = set()
    rows_examined = retained_rows = 0

    for row in feature_rows:
        rows_examined += 1
        for field in EXPECTED_CATEGORIES:
            categories[field].add(row.get(field, ""))
        pair = (row.get("sample", ""), row.get("virus", ""))
        if pair not in primary_pairs:
            continue
        if row.get("mapping_mode") != "exact" or row.get("virus_assignment") != "assigned":
            continue
        strand = row.get("strand", "")
        if strand not in STRANDS:
            continue
        retained_rows += 1
        try:
            declared_length = int(row["length"])
            count = float(row["count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise Stage03Error(f"invalid numeric read-level row {rows_examined}: {exc}") from exc
        sequence = row.get("sequence", "")
        if declared_length != len(sequence):
            raise Stage03Error(
                f"declared/actual length mismatch at read-level row {rows_examined}: "
                f"{declared_length}!={len(sequence)}"
            )
        if count < 0 or not math.isfinite(count):
            raise Stage03Error(f"invalid count at read-level row {rows_examined}: {count}")
        if set(sequence) - set("ACGT"):
            unexpected_sequences.add(sequence)
        if config.passenger_length_min <= declared_length <= config.passenger_length_max:
            sequence_abundance[(pair[0], pair[1], strand, declared_length)][sequence] += count

    runs: list[dict[str, object]] = []
    focal_manifest: list[dict[str, object]] = []
    passenger_manifest: list[dict[str, object]] = []
    focal_by_run: dict[str, list[dict[str, object]]] = {}
    passenger_by_run: dict[str, list[dict[str, object]]] = {}

    for sample, unit in sorted(primary_pairs):
        pair_meta = metadata[(sample, unit)]
        for focal_length in config.focal_lengths:
            for focal_strand in STRANDS:
                run_id = run_identifier(sample, unit, focal_length, focal_strand)
                focals = []
                focal_sequences = sequence_abundance.get(
                    (sample, unit, focal_strand, focal_length), Counter()
                )
                for sequence, abundance in sorted(focal_sequences.items()):
                    row = {
                        "focal_id": stable_identifier(
                            "F", sample, unit, focal_length, focal_strand, sequence
                        ),
                        "sequence": sequence,
                        "sample": sample,
                        "analysis_unit": unit,
                        "biological_virus": pair_meta.get("biological_virus", ""),
                        "focal_length": focal_length,
                        "focal_strand": focal_strand,
                        "focal_abundance": abundance,
                        "run_id": run_id,
                    }
                    focals.append(row)
                    focal_manifest.append(row)
                passengers = []
                passenger_strand = opposite_strand(focal_strand)
                for length in range(
                    config.passenger_length_min, config.passenger_length_max + 1
                ):
                    passenger_sequences = sequence_abundance.get(
                        (sample, unit, passenger_strand, length), Counter()
                    )
                    for sequence, abundance in sorted(passenger_sequences.items()):
                        row = {
                            "passenger_id": stable_identifier("P", run_id, sequence),
                            "sequence": sequence,
                            "sample": sample,
                            "analysis_unit": unit,
                            "biological_virus": pair_meta.get("biological_virus", ""),
                            "focal_length": focal_length,
                            "focal_strand": focal_strand,
                            "passenger_strand": passenger_strand,
                            "passenger_length": length,
                            "passenger_abundance": abundance,
                            "run_id": run_id,
                        }
                        passengers.append(row)
                        passenger_manifest.append(row)
                focal_by_run[run_id] = focals
                passenger_by_run[run_id] = passengers
                runs.append({
                    "run_id": run_id,
                    "sample": sample,
                    "analysis_unit": unit,
                    "biological_virus": pair_meta.get("biological_virus", ""),
                    "polarity": pair_meta.get("polarity", ""),
                    "focal_length": focal_length,
                    "focal_strand": focal_strand,
                    "passenger_strand": passenger_strand,
                    "n_focal_references": len(focals),
                    "total_focal_abundance": sum(float(x["focal_abundance"]) for x in focals),
                    "n_passenger_candidates": len(passengers),
                })
    return {
        "runs": runs,
        "focal_manifest": focal_manifest,
        "passenger_manifest": passenger_manifest,
        "focal_by_run": focal_by_run,
        "passenger_by_run": passenger_by_run,
        "rows_examined": rows_examined,
        "retained_rows": retained_rows,
        "categories": categories,
        "unexpected_sequences": unexpected_sequences,
        "primary_pairs": primary_pairs,
    }


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "NOT_INSTALLED"


def command_version(command: list[str]) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    return (result.stdout or result.stderr).splitlines()[0].strip()


def official_distances(bam, alignment) -> tuple[int, int]:
    from stepRNA.commands import left_overhang, right_overhang

    positions = alignment.get_reference_positions(full_length=True)
    left, _ = left_overhang(bam, alignment, positions)
    right, _ = right_overhang(bam, alignment, positions)
    return int(left), int(right)


def run_steprna(reference: Path, reads: Path, name: str, directory: Path) -> subprocess.CompletedProcess[str]:
    directory.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        [
            "stepRNA", "--reference", str(reference), "--reads", str(reads),
            "--name", name, "--directory", str(directory),
        ],
        text=True,
        capture_output=True,
    )


def synthetic_records() -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    from Bio.Seq import Seq

    refs = {
        "FJ": "ACGTTGCACTGATCGTACGATGC",
        "FB": "GATCCGTAGTCAGGCTAACGTCA",
        "FO": "TGCAGATCTAGCGTACCTGATCA",
    }
    reads = {
        "QJ": str(Seq("TT" + refs["FJ"][:21]).reverse_complement()),
        "QB": str(Seq(refs["FB"]).reverse_complement()),
        "QO": str(Seq(refs["FO"][2:]).reverse_complement()),
        "QD": "AAAAAAAAAAAAAAA",
    }
    return list(refs.items()), list(reads.items())


def run_preflight(output_root: Path, required_version: str = "1.0.6") -> dict[str, object]:
    import pysam

    preflight = output_root / "provenance" / "preflight"
    if preflight.exists():
        shutil.rmtree(preflight)
    preflight.mkdir(parents=True)
    reference = preflight / "synthetic_fileA.fa"
    reads = preflight / "synthetic_fileB.fa"
    write_fasta(reference, synthetic_records()[0])
    write_fasta(reads, synthetic_records()[1])

    steprna_version = command_version(["stepRNA", "--version"]).removeprefix("stepRNA v")
    bowtie2_version = command_version(["bowtie2", "--version"])
    versions = {
        "stepRNA": steprna_version,
        "Bowtie2": bowtie2_version,
        "Python": sys.version.split()[0],
        "pysam": package_version("pysam"),
        "Biopython": package_version("biopython"),
        "NumPy": package_version("numpy"),
        "alive-progress": package_version("alive-progress"),
    }
    if steprna_version != required_version:
        raise Stage03Error(f"stepRNA version {steprna_version}; required {required_version}")

    official_dir = preflight / "official"
    result = run_steprna(reference, reads, "synthetic", official_dir)
    (preflight / "stdout.txt").write_text(result.stdout)
    (preflight / "stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise Stage03Error(f"official synthetic stepRNA failed with exit code {result.returncode}")
    expected = [official_dir / f"synthetic{suffix}" for suffix in REQUIRED_OFFICIAL_SUFFIXES]
    expected += [
        official_dir / "synthetic.sorted.bam",
        official_dir / "synthetic_AlignmentFiles" / "synthetic_passed.bam",
    ]
    missing = [str(x) for x in expected if not x.exists()]
    if missing:
        raise Stage03Error(f"official synthetic outputs missing: {','.join(missing)}")
    indexes = list(preflight.glob("synthetic_fileA.*.bt2"))
    if len(indexes) < 6:
        raise Stage03Error("official synthetic Bowtie2 index was not produced")

    recovered: dict[tuple[str, str], tuple[int, int]] = {}
    passed = official_dir / "synthetic_AlignmentFiles" / "synthetic_passed.bam"
    with pysam.AlignmentFile(passed, "rb") as bam:
        for alignment in bam:
            recovered[(alignment.reference_name, alignment.query_name)] = official_distances(bam, alignment)
    required = {
        ("FJ", "QJ"): (2, -2),
        ("FB", "QB"): (0, 0),
        ("FO", "QO"): (-2, 0),
    }
    for key, geometry in required.items():
        if recovered.get(key) != geometry:
            raise Stage03Error(f"synthetic geometry {key}: {recovered.get(key)} != {geometry}")

    checks = [
        {"check": "official_executable_launch", "status": "PASS", "details": steprna_version},
        {"check": "bowtie2_index_alignment", "status": "PASS", "details": bowtie2_version},
        {"check": "official_summary_files", "status": "PASS", "details": str(len(expected))},
        {"check": "classified_bam_parse", "status": "PASS", "details": str(len(recovered))},
        {"check": "negative_overhang_sign", "status": "PASS", "details": "FO/QO 5p=-2"},
        {"check": "zero_blunt_sign", "status": "PASS", "details": "FB/QB=(0,0)"},
        {"check": "positive_underhang_sign", "status": "PASS", "details": "FJ/QJ 5p=+2"},
        {"check": "known_joint_geometry", "status": "PASS", "details": "FJ/QJ=(+2,-2)"},
    ]
    write_table(
        output_root / "provenance" / "software_versions.tsv",
        [{"software": k, "version": v} for k, v in versions.items()],
        ["software", "version"],
    )
    write_table(
        output_root / "provenance" / "preflight.tsv",
        checks,
        ["check", "status", "details"],
    )
    return {"versions": versions, "checks": checks, "geometries": recovered}


def parse_number(value: str) -> float | None:
    return None if value in {"", "NA", "nan", "NaN"} else float(value)


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))
    except (OSError, csv.Error) as exc:
        raise Stage03Error(f"cannot parse official file {path}: {exc}") from exc


def summarise_duplexes(
    run: dict[str, object],
    focal_rows: list[dict[str, object]],
    duplexes: list[dict[str, object]],
    config: Stage03Config,
) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]]]:
    focal_by_id = {str(row["focal_id"]): row for row in focal_rows}
    recovered_counts = Counter(str(x["focal_id"]) for x in duplexes)
    recovered_ids = set(recovered_counts)
    unknown = recovered_ids - set(focal_by_id)
    if unknown:
        raise Stage03Error(f"duplex focal IDs absent from manifest: {','.join(sorted(unknown))}")
    joint = [
        x for x in duplexes
        if x["steprna_5p_distance"] == config.joint_5p_distance
        and x["steprna_3p_distance"] == config.joint_3p_distance
    ]
    joint_ids = {str(x["focal_id"]) for x in joint}
    total_abundance = sum(float(x["focal_abundance"]) for x in focal_rows)
    recovered_abundance = sum(float(focal_by_id[x]["focal_abundance"]) for x in recovered_ids)
    joint_abundance = sum(float(focal_by_id[x]["focal_abundance"]) for x in joint_ids)
    common = {k: run[k] for k in (
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand"
    )}
    recovery = {
        **common,
        "n_focal_references": len(focal_rows),
        "n_recovered_focal_references": len(recovered_ids),
        "passenger_recovery_fraction_unique": safe_fraction(len(recovered_ids), len(focal_rows)),
        "total_focal_abundance": total_abundance,
        "recovered_focal_abundance": recovered_abundance,
        "passenger_recovery_fraction_abundance": safe_fraction(recovered_abundance, total_abundance),
        "run_id": run["run_id"],
    }
    geometry = {
        **common,
        "steprna_5p_distance": config.joint_5p_distance,
        "steprna_3p_distance": config.joint_3p_distance,
        "n_recovered_duplexes": len(duplexes),
        "n_joint_geometry_duplexes": len(joint),
        "varroa_2nt_joint_duplex_fraction": safe_fraction(len(joint), len(duplexes)),
        "n_focal_references": len(focal_rows),
        "n_recovered_focal_references": len(recovered_ids),
        "n_focal_references_supporting_joint_geometry": len(joint_ids),
        "varroa_2nt_reference_fraction_all": safe_fraction(len(joint_ids), len(focal_rows)),
        "varroa_2nt_reference_fraction_recovered": safe_fraction(len(joint_ids), len(recovered_ids)),
        "total_focal_abundance": total_abundance,
        "recovered_focal_abundance": recovered_abundance,
        "joint_supporting_focal_abundance": joint_abundance,
        "varroa_2nt_reference_fraction_abundance_all": safe_fraction(joint_abundance, total_abundance),
        "varroa_2nt_reference_fraction_abundance_recovered": safe_fraction(joint_abundance, recovered_abundance),
        "run_id": run["run_id"],
    }
    joint_counts = Counter(str(x["focal_id"]) for x in joint)
    joint_references = []
    for focal_id in sorted(joint_ids):
        focal = focal_by_id[focal_id]
        joint_references.append({
            **common,
            "focal_id": focal_id,
            "sequence": focal["sequence"],
            "focal_abundance": focal["focal_abundance"],
            "n_recovered_duplexes_for_reference": recovered_counts[focal_id],
            "n_joint_geometry_duplexes_for_reference": joint_counts[focal_id],
            "run_id": run["run_id"],
        })
    return recovery, geometry, joint_references


def parse_official_run(
    run: dict[str, object],
    raw_dir: Path,
    focal_rows: list[dict[str, object]],
    passenger_rows: list[dict[str, object]],
    config: Stage03Config,
) -> dict[str, object]:
    import pysam

    run_id = str(run["run_id"])
    paths = {suffix: raw_dir / f"{run_id}{suffix}" for suffix in REQUIRED_OFFICIAL_SUFFIXES}
    passed_bam = raw_dir / f"{run_id}_AlignmentFiles" / f"{run_id}_passed.bam"
    missing = [str(path) for path in [*paths.values(), passed_bam] if not path.exists()]
    if missing:
        raise Stage03Error(f"missing official files for {run_id}: {','.join(missing)}")

    focal_by_id = {str(row["focal_id"]): row for row in focal_rows}
    passenger_by_id = {str(row["passenger_id"]): row for row in passenger_rows}
    duplexes: list[dict[str, object]] = []
    try:
        with pysam.AlignmentFile(passed_bam, "rb") as bam:
            for alignment in bam:
                focal_id = alignment.reference_name
                passenger_id = alignment.query_name
                if focal_id not in focal_by_id or passenger_id not in passenger_by_id:
                    raise Stage03Error(
                        f"manifest/parser identifier mismatch in {run_id}: {focal_id}/{passenger_id}"
                    )
                d5, d3 = official_distances(bam, alignment)
                duplexes.append({
                    "focal_id": focal_id,
                    "passenger_id": passenger_id,
                    "steprna_5p_distance": d5,
                    "steprna_3p_distance": d3,
                    "passenger_length": alignment.query_length,
                })
    except (OSError, ValueError) as exc:
        raise Stage03Error(f"cannot parse classified BAM for {run_id}: {exc}") from exc

    overhang = read_csv(paths["_overhang.csv"])
    unique = {int(row["Overhang"]): row for row in read_csv(paths["_unique_overhang.csv"])}
    spectrum = []
    for row in overhang:
        try:
            distance = int(row["Overhang"])
            unique_row = unique[distance]
        except (KeyError, TypeError, ValueError) as exc:
            raise Stage03Error(f"malformed official overhang row in {run_id}: {row}") from exc
        for end, prefix in (("5p", "5prime"), ("3p", "3prime")):
            spectrum.append({
                **{k: run[k] for k in (
                    "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand"
                )},
                "end": end,
                "signed_distance": distance,
                "official_duplex_count": int(row[prefix]),
                "official_unique_reference_count": int(unique_row[prefix]),
                "official_duplex_log_ratio": parse_number(row[f"{prefix}_Logodds"]),
                "official_duplex_wald_z": parse_number(row[f"{prefix}_Zscore"]),
                "official_unique_reference_log_ratio": parse_number(unique_row[f"{prefix}_Logodds"]),
                "official_unique_reference_wald_z": parse_number(unique_row[f"{prefix}_Zscore"]),
                "run_id": run_id,
            })

    passenger_length_rows = []
    bam_lengths = Counter(int(x["passenger_length"]) for x in duplexes)
    bam_length_refs: dict[int, set[str]] = defaultdict(set)
    for x in duplexes:
        bam_length_refs[int(x["passenger_length"])].add(str(x["focal_id"]))
    official_lengths = {}
    for row in read_csv(paths["_passenger_length.csv"]):
        try:
            length, count = int(row["passenger_length"]), int(row["passenger_count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise Stage03Error(f"malformed passenger-length row in {run_id}: {row}") from exc
        official_lengths[length] = count
        passenger_length_rows.append({
            **{k: run[k] for k in (
                "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand"
            )},
            "passenger_length": length,
            "official_duplex_count": count,
            "official_unique_reference_count": len(bam_length_refs[length]),
            "run_id": run_id,
        })
    if official_lengths != dict(bam_lengths):
        raise Stage03Error(f"passenger-length/BAM inconsistency in {run_id}")

    passenger_numbers = {}
    for row in read_csv(paths["_passenger_number.csv"]):
        try:
            passenger_numbers[row["siRNA_reference"]] = int(row["number"])
        except (KeyError, TypeError, ValueError) as exc:
            raise Stage03Error(f"malformed passenger-number row in {run_id}: {row}") from exc
    bam_reference_counts = Counter(str(x["focal_id"]) for x in duplexes)
    expected_reference_counts = {focal_id: bam_reference_counts[focal_id] for focal_id in focal_by_id}
    if passenger_numbers != expected_reference_counts:
        raise Stage03Error(f"passenger-number/BAM inconsistency in {run_id}")

    for end, field in (("5p", "steprna_5p_distance"), ("3p", "steprna_3p_distance")):
        official_counts = {
            int(row["signed_distance"]): int(row["official_duplex_count"])
            for row in spectrum if row["end"] == end
        }
        bam_counts = Counter(int(x[field]) for x in duplexes)
        if any(official_counts.get(key, 0) != value for key, value in bam_counts.items()):
            raise Stage03Error(f"{end} overhang/BAM inconsistency in {run_id}")
        if sum(official_counts.values()) != len(duplexes):
            raise Stage03Error(f"{end} total/BAM inconsistency in {run_id}")

    recovery, geometry, joint_references = summarise_duplexes(
        run, focal_rows, duplexes, config
    )
    return {
        "spectrum": spectrum,
        "passenger_lengths": passenger_length_rows,
        "recovery": recovery,
        "geometry": geometry,
        "joint_references": joint_references,
        "duplexes": duplexes,
    }


def zero_run_summaries(
    run: dict[str, object], focal_rows: list[dict[str, object]]
) -> tuple[dict[str, object], dict[str, object]]:
    total_abundance = sum(float(x["focal_abundance"]) for x in focal_rows)
    common = {k: run[k] for k in (
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand"
    )}
    recovery = {
        **common,
        "n_focal_references": len(focal_rows),
        "n_recovered_focal_references": 0,
        "passenger_recovery_fraction_unique": safe_fraction(0, len(focal_rows)),
        "total_focal_abundance": total_abundance,
        "recovered_focal_abundance": 0,
        "passenger_recovery_fraction_abundance": safe_fraction(0, total_abundance),
        "run_id": run["run_id"],
    }
    geometry = {
        **common,
        "steprna_5p_distance": 2,
        "steprna_3p_distance": -2,
        "n_recovered_duplexes": 0,
        "n_joint_geometry_duplexes": 0,
        "varroa_2nt_joint_duplex_fraction": None,
        "n_focal_references": len(focal_rows),
        "n_recovered_focal_references": 0,
        "n_focal_references_supporting_joint_geometry": 0,
        "varroa_2nt_reference_fraction_all": safe_fraction(0, len(focal_rows)),
        "varroa_2nt_reference_fraction_recovered": None,
        "total_focal_abundance": total_abundance,
        "recovered_focal_abundance": 0,
        "joint_supporting_focal_abundance": 0,
        "varroa_2nt_reference_fraction_abundance_all": safe_fraction(0, total_abundance),
        "varroa_2nt_reference_fraction_abundance_recovered": None,
        "run_id": run["run_id"],
    }
    return recovery, geometry


FIELDS = {
    "input_manifest": [
        "run_id", "sample", "analysis_unit", "biological_virus", "polarity",
        "focal_length", "focal_strand", "passenger_strand", "file_a", "file_b",
        "n_focal_references", "total_focal_abundance", "n_passenger_candidates",
    ],
    "focal_manifest": [
        "focal_id", "sequence", "sample", "analysis_unit", "biological_virus",
        "focal_length", "focal_strand", "focal_abundance", "run_id",
    ],
    "passenger_manifest": [
        "passenger_id", "sequence", "sample", "analysis_unit", "biological_virus",
        "focal_length", "focal_strand", "passenger_strand", "passenger_length",
        "passenger_abundance", "run_id",
    ],
    "run_manifest": [
        "run_id", "sample", "analysis_unit", "biological_virus", "focal_length",
        "focal_strand", "passenger_strand", "n_focal_references",
        "n_passenger_candidates", "status", "exit_code", "runtime_seconds",
        "raw_output_directory", "message",
    ],
    "recovery": [
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand",
        "n_focal_references", "n_recovered_focal_references",
        "passenger_recovery_fraction_unique", "total_focal_abundance",
        "recovered_focal_abundance", "passenger_recovery_fraction_abundance", "run_id",
    ],
    "spectrum": [
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand",
        "end", "signed_distance", "official_duplex_count",
        "official_unique_reference_count", "official_duplex_log_ratio",
        "official_duplex_wald_z", "official_unique_reference_log_ratio",
        "official_unique_reference_wald_z", "run_id",
    ],
    "passenger_lengths": [
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand",
        "passenger_length", "official_duplex_count",
        "official_unique_reference_count", "run_id",
    ],
    "geometry": [
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand",
        "steprna_5p_distance", "steprna_3p_distance", "n_recovered_duplexes",
        "n_joint_geometry_duplexes", "varroa_2nt_joint_duplex_fraction",
        "n_focal_references", "n_recovered_focal_references",
        "n_focal_references_supporting_joint_geometry",
        "varroa_2nt_reference_fraction_all", "varroa_2nt_reference_fraction_recovered",
        "total_focal_abundance", "recovered_focal_abundance",
        "joint_supporting_focal_abundance",
        "varroa_2nt_reference_fraction_abundance_all",
        "varroa_2nt_reference_fraction_abundance_recovered", "run_id",
    ],
    "joint_references": [
        "sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand",
        "focal_id", "sequence", "focal_abundance",
        "n_recovered_duplexes_for_reference",
        "n_joint_geometry_duplexes_for_reference", "run_id",
    ],
    "qc": ["metric", "status", "value", "details"],
}


def run_stage03(legacy_core: Path, config_path: Path, output_root: Path) -> tuple[float, bool]:
    started = time.monotonic()
    config = load_config(config_path)
    preflight = run_preflight(output_root, config.required_steprna_version)

    eligibility = read_eligibility(legacy_core / "results/descriptive/eligibility.tsv")
    samples = {row["sample"] for row in eligibility}
    inputs = collapse_inputs(eligibility, iter_feature_rows(legacy_core, samples), config)
    runs = list(inputs["runs"])
    focal_by_run = inputs["focal_by_run"]
    passenger_by_run = inputs["passenger_by_run"]
    input_dir = output_root / "inputs" / "fasta"

    input_manifest = []
    for run in runs:
        run_id = str(run["run_id"])
        file_a = input_dir / f"{run_id}.fileA.fa"
        file_b = input_dir / f"{run_id}.fileB.fa"
        write_fasta(file_a, ((str(x["focal_id"]), str(x["sequence"])) for x in focal_by_run[run_id]))
        write_fasta(file_b, ((str(x["passenger_id"]), str(x["sequence"])) for x in passenger_by_run[run_id]))
        input_manifest.append({**run, "file_a": str(file_a), "file_b": str(file_b)})

    write_table(output_root / "inputs" / "input_manifest.tsv", input_manifest, FIELDS["input_manifest"])
    write_table(
        output_root / "inputs" / "focal_reference_manifest.tsv.gz",
        list(inputs["focal_manifest"]), FIELDS["focal_manifest"], gzip_output=True,
    )
    write_table(
        output_root / "inputs" / "passenger_manifest.tsv.gz",
        list(inputs["passenger_manifest"]), FIELDS["passenger_manifest"], gzip_output=True,
    )

    run_manifest = []
    recoveries = []
    spectra = []
    passenger_lengths = []
    geometries = []
    joint_references = []
    malformed = missing_files = id_inconsistencies = geometry_inconsistencies = 0

    for run, input_row in zip(runs, input_manifest):
        run_started = time.monotonic()
        run_id = str(run["run_id"])
        focals = focal_by_run[run_id]
        passengers = passenger_by_run[run_id]
        raw_dir = output_root / "raw" / run_id
        status = "success"
        exit_code: int | None = 0
        message = ""
        if not focals:
            status = "zero_focal"
            recovery, geometry = zero_run_summaries(run, focals)
            recoveries.append(recovery)
            geometries.append(geometry)
        elif not passengers:
            status = "zero_passenger_pool"
            recovery, geometry = zero_run_summaries(run, focals)
            recoveries.append(recovery)
            geometries.append(geometry)
        else:
            if raw_dir.exists():
                shutil.rmtree(raw_dir)
            raw_dir.mkdir(parents=True)
            official_a = raw_dir / "fileA.fa"
            official_b = raw_dir / "fileB.fa"
            shutil.copy2(input_row["file_a"], official_a)
            shutil.copy2(input_row["file_b"], official_b)
            result = run_steprna(official_a, official_b, run_id, raw_dir)
            (raw_dir / "official_stdout.txt").write_text(result.stdout)
            (raw_dir / "official_stderr.txt").write_text(result.stderr)
            exit_code = result.returncode
            if result.returncode:
                status = "failed"
                message = f"official stepRNA exit code {result.returncode}"
            else:
                try:
                    parsed = parse_official_run(run, raw_dir, focals, passengers, config)
                    recoveries.append(parsed["recovery"])
                    geometries.append(parsed["geometry"])
                    spectra.extend(parsed["spectrum"])
                    passenger_lengths.extend(parsed["passenger_lengths"])
                    joint_references.extend(parsed["joint_references"])
                except Stage03Error as exc:
                    status = "failed"
                    message = str(exc)
                    lower = message.lower()
                    missing_files += int("missing official" in lower)
                    id_inconsistencies += int("identifier mismatch" in lower)
                    geometry_inconsistencies += int("inconsistency" in lower)
                    malformed += int("malformed" in lower or "cannot parse" in lower)
        run_manifest.append({
            **{k: run[k] for k in (
                "run_id", "sample", "analysis_unit", "biological_virus", "focal_length",
                "focal_strand", "passenger_strand", "n_focal_references",
                "n_passenger_candidates",
            )},
            "status": status,
            "exit_code": exit_code,
            "runtime_seconds": time.monotonic() - run_started,
            "raw_output_directory": str(raw_dir) if status == "success" else "NA",
            "message": message,
        })

    write_table(output_root / "provenance" / "run_manifest.tsv", run_manifest, FIELDS["run_manifest"])
    write_table(output_root / "parsed" / "passenger_recovery_by_pair.tsv", recoveries, FIELDS["recovery"])
    write_table(output_root / "parsed" / "overhang_spectrum_by_pair.tsv", spectra, FIELDS["spectrum"])
    write_table(output_root / "parsed" / "passenger_length_by_pair.tsv", passenger_lengths, FIELDS["passenger_lengths"])
    write_table(output_root / "parsed" / "joint_geometry_by_pair.tsv", geometries, FIELDS["geometry"])
    write_table(
        output_root / "parsed" / "joint_geometry_references.tsv.gz",
        joint_references, FIELDS["joint_references"], gzip_output=True,
    )

    qc: list[dict[str, object]] = []
    def q(metric: str, value: object, status: str = "PASS", details: str = "") -> None:
        qc.append({"metric": metric, "status": status, "value": value, "details": details})

    primary_pairs = inputs["primary_pairs"]
    q("primary_eligible_samples", len({x[0] for x in primary_pairs}))
    q("primary_eligible_sample_virus_units", len(primary_pairs))
    q("maximum_possible_focal_populations", len(primary_pairs) * 4)
    q("non_zero_focal_populations", sum(bool(focal_by_run[str(x["run_id"])]) for x in runs))
    q("zero_focal_populations", sum(not focal_by_run[str(x["run_id"])] for x in runs))
    for length in config.focal_lengths:
        for strand in STRANDS:
            relevant = [
                x for x in inputs["focal_manifest"]
                if x["focal_length"] == length and x["focal_strand"] == strand
            ]
            q("distinct_focal_sequences", len(relevant), "INFO", f"length={length}; strand={strand}")
            q("represented_focal_abundance", sum(float(x["focal_abundance"]) for x in relevant), "INFO", f"length={length}; strand={strand}")
    q("distinct_passenger_sequences_across_runs", len(inputs["passenger_manifest"]))
    passenger_input_lengths = [int(x["passenger_length"]) for x in inputs["passenger_manifest"]]
    q("passenger_input_length_min", min(passenger_input_lengths, default="NA"))
    q("passenger_input_length_max", max(passenger_input_lengths, default="NA"))
    q("read_level_rows_examined", inputs["rows_examined"])
    q("canonical_exact_assigned_rows_retained", inputs["retained_rows"])
    for field, allowed in EXPECTED_CATEGORIES.items():
        unexpected = sorted(inputs["categories"][field] - allowed)
        q(f"unexpected_{field}_categories", len(unexpected), "WARN" if unexpected else "PASS", ",".join(unexpected))
    q("unexpected_sequence_alphabet_values", len(inputs["unexpected_sequences"]), "WARN" if inputs["unexpected_sequences"] else "PASS")
    q("official_steprna_version", preflight["versions"]["stepRNA"])
    q("bowtie2_version", preflight["versions"]["Bowtie2"])
    q("successful_biological_runs", sum(x["status"] == "success" for x in run_manifest))
    q("failed_biological_runs", sum(x["status"] == "failed" for x in run_manifest), "FAIL" if any(x["status"] == "failed" for x in run_manifest) else "PASS")
    q("zero_passenger_pool_runs", sum(x["status"] == "zero_passenger_pool" for x in run_manifest))
    zero_recovered = sum(
        x["status"] == "success" and next(
            (g["n_recovered_duplexes"] for g in geometries if g["run_id"] == x["run_id"]), 0
        ) == 0 for x in run_manifest
    )
    q("successful_runs_zero_recovered_passengers", zero_recovered)
    q("number_of_focal_references", sum(int(x["n_focal_references"]) for x in recoveries))
    q("focal_references_with_passengers", sum(int(x["n_recovered_focal_references"]) for x in recoveries))
    q("total_reconstructed_duplex_relationships", sum(int(x["n_recovered_duplexes"]) for x in geometries))
    distances = [int(x["signed_distance"]) for x in spectra if int(x["official_duplex_count"]) > 0]
    q("signed_distance_min", min(distances, default="NA"))
    q("signed_distance_max", max(distances, default="NA"))
    parsed_lengths = [int(x["passenger_length"]) for x in passenger_lengths if int(x["official_duplex_count"]) > 0]
    q("passenger_length_min", min(parsed_lengths, default="NA"))
    q("passenger_length_max", max(parsed_lengths, default="NA"))
    q("malformed_unparseable_official_rows", malformed, "FAIL" if malformed else "PASS")
    q("runs_missing_official_files", missing_files, "FAIL" if missing_files else "PASS")
    q("focal_id_manifest_parser_inconsistencies", id_inconsistencies, "FAIL" if id_inconsistencies else "PASS")
    q("joint_geometry_alignment_inconsistencies", geometry_inconsistencies, "FAIL" if geometry_inconsistencies else "PASS")
    na_counts = sum(
        value is None
        for row in geometries
        for key, value in row.items()
        if key.startswith("varroa_2nt_")
    )
    q("na_denominator_metric_values", na_counts, "INFO")
    write_table(output_root / "qc" / "stage03_accounting.tsv", qc, FIELDS["qc"])
    failed = any(x["status"] == "FAIL" for x in qc)
    return time.monotonic() - started, failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args(argv)
    started = time.monotonic()
    try:
        config = load_config(args.config.resolve())
        if args.preflight_only:
            run_preflight(args.output_root.resolve(), config.required_steprna_version)
            print(f"Stage 03 preflight completed in {time.monotonic() - started:.3f} seconds", file=sys.stderr)
            return 0
        if args.legacy_core is None:
            parser.error("--legacy-core is required unless --preflight-only is used")
        elapsed, failed = run_stage03(
            args.legacy_core.resolve(), args.config.resolve(), args.output_root.resolve()
        )
        print(f"Stage 03 completed in {elapsed:.3f} seconds", file=sys.stderr)
        return 1 if failed else 0
    except Exception as exc:
        output_root = args.output_root.resolve()
        failure = [{"metric": "stage03_execution", "status": "FAIL", "value": 1, "details": str(exc)}]
        write_table(output_root / "qc" / "stage03_accounting.tsv", failure, FIELDS["qc"])
        print(f"Stage 03 failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
