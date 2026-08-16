#!/usr/bin/env python3
"""Canonical Stage 09A training-universe construction and orchestration.

The Python layer owns frozen-input accounting, exact sequence predictors, and
deterministic data preparation.  The approved R/glmnet layer owns penalized
model fitting and validation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import statistics
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


LENGTHS = (23, 24)
DNA = frozenset("ACGT")
TERMINAL_POSITIONS = ("5p1", "5p2", "3p2", "3p1")
BASE_FEATURE_NAMES = tuple(
    f"guide_{position}_{base}"
    for position in TERMINAL_POSITIONS
    for base in ("C", "G", "U")
) + ("guide_A3p3", "guide_GC_3p5_10", "guide_W17", "guide_R10")
EXPECTED_ACCOUNTING = {
    "primary_samples": 20,
    "sample_virus_units": 54,
    "sample_virus_length_groups": 108,
    "opportunities_23nt": 411079,
    "opportunities_24nt": 408148,
    "opportunities_total": 819227,
    "represented_23nt": 121592,
    "represented_24nt": 175564,
    "represented_total": 297156,
    "supported_abundance": 3445943.0,
}
ALPHA_GRID = (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0)
L1_RATIO_GRID = (0.05, 0.25, 0.50, 0.75, 0.95, 1.0)
STRUCTURES = ("A_shared", "B_shared_plus_length_interactions", "C_separate_23_24")
STAGE08_FORBIDDEN = {
    "target_whole_p_unpaired",
    "target_seed_g2_8_p_unpaired",
    "guide_5p_terminal_dg_4bp",
    "passenger_5p_terminal_dg_4bp",
    "asymmetry_ddg_4bp",
    "guide_5p_terminal_dg_5bp",
    "passenger_5p_terminal_dg_5bp",
    "asymmetry_ddg_5bp",
    "guide_self_fold_mfe_kcal_mol",
    "guide_self_fold_structure",
}
COMPLEMENT = str.maketrans("ACGTU", "TGCAA")


@dataclass(frozen=True)
class Scaling:
    means: dict[str, float]
    sds: dict[str, float]
    retained: tuple[str, ...]
    omitted: tuple[str, ...]


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def reverse_complement_dna(sequence: str) -> str:
    normalized = sequence.upper().replace("U", "T")
    if not normalized or any(base not in DNA for base in normalized):
        raise ValueError("reverse complement requires a non-empty A/C/G/T(U) sequence")
    return normalized.translate(COMPLEMENT)[::-1]


def reverse_complement_rna(sequence: str) -> str:
    return reverse_complement_dna(sequence).replace("T", "U")


def guide_predictors(sequence: str) -> dict[str, object]:
    """Return the eight physical-guide predictors in v0.19.1 orientation."""
    guide = sequence.upper().replace("T", "U")
    if len(guide) not in LENGTHS or any(base not in "ACGU" for base in guide):
        raise ValueError("Stage 09A guide must be an A/C/G/U 23- or 24-mer")
    gc_slice = guide[-10:-4]
    if len(gc_slice) != 6:
        raise AssertionError("GC_3p5_10 slice must contain exactly six bases")
    return {
        "guide_5p1_nt": guide[0],
        "guide_5p2_nt": guide[1],
        "guide_3p2_nt": guide[-2],
        "guide_3p1_nt": guide[-1],
        "guide_A3p3": int(guide[-3] == "A"),
        "guide_GC_3p5_10": sum(base in "GC" for base in gc_slice) / 6.0,
        "guide_W17": int(guide[16] in "AU"),
        "guide_R10": int(guide[9] in "AG"),
    }


def encode_predictors(predictors: Mapping[str, object]) -> dict[str, float]:
    encoded: dict[str, float] = {}
    source = {
        "5p1": str(predictors["guide_5p1_nt"]),
        "5p2": str(predictors["guide_5p2_nt"]),
        "3p2": str(predictors["guide_3p2_nt"]),
        "3p1": str(predictors["guide_3p1_nt"]),
    }
    for position, nucleotide in source.items():
        if nucleotide not in "ACGU":
            raise ValueError(f"unexpected terminal nucleotide {nucleotide!r}")
        for category in "CGU":
            encoded[f"guide_{position}_{category}"] = float(nucleotide == category)
    for name in ("guide_A3p3", "guide_GC_3p5_10", "guide_W17", "guide_R10"):
        encoded[name] = float(predictors[name])
    if tuple(encoded) != BASE_FEATURE_NAMES:
        raise AssertionError("encoded Stage 09A feature order changed")
    return encoded


def parse_fasta(path: Path) -> list[str]:
    records: list[str] = []
    sequence: list[str] = []
    seen = False
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen:
                    records.append("".join(sequence).upper().replace("U", "T"))
                seen, sequence = True, []
            else:
                if not seen:
                    raise ValueError(f"sequence before FASTA header: {path}:{line_number}")
                sequence.append(line)
    if seen:
        records.append("".join(sequence).upper().replace("U", "T"))
    if not records or any(not record for record in records):
        raise ValueError(f"missing or empty FASTA record: {path}")
    return records


def supported_antisense_sequences(records: Iterable[str], length: int) -> set[str]:
    """Enumerate unique depth-supported antisense sequence opportunities."""
    if length not in LENGTHS:
        raise ValueError("canonical Stage 09A supports only 23 and 24 nt")
    output: set[str] = set()
    for record in records:
        for start in range(max(0, len(record) - length + 1)):
            window = record[start : start + length]
            if all(base in DNA for base in window):
                output.add(reverse_complement_dna(window))
    return output


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_primary_pairs(legacy_core: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = legacy_core / "results/descriptive/eligibility.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"missing frozen eligibility table: {path}")
    rows = _read_tsv(path)
    required = {"sample", "analysis_unit", "biological_virus", "primary_eligible"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("frozen eligibility table lacks required Stage 09A columns")
    pairs = {
        (row["sample"], row["analysis_unit"]): row
        for row in rows
        if is_true(row["primary_eligible"])
    }
    if len(pairs) != len(set(pairs)):
        raise ValueError("duplicate primary-eligible sample-virus pair")
    return pairs


def aggregate_observed_antisense(
    legacy_core: Path, pairs: Mapping[tuple[str, str], Mapping[str, str]]
) -> tuple[dict[tuple[str, str, int], Counter[str]], dict[str, int]]:
    """Sum canonical abundance by physical antisense sequence."""
    output: dict[tuple[str, str, int], Counter[str]] = defaultdict(Counter)
    accounting = Counter()
    samples = sorted({sample for sample, _unit in pairs})
    expected_columns = {
        "sample", "mapping_mode", "virus", "virus_assignment", "strand",
        "sequence", "length", "count",
    }
    for sample in samples:
        path = legacy_core / "tables" / sample / f"{sample}.read_level_features.tsv.gz"
        if not path.is_file():
            raise FileNotFoundError(f"missing frozen read-level table: {path}")
        try:
            with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None or not expected_columns.issubset(reader.fieldnames):
                    raise ValueError(f"invalid Stage 09A read-level schema: {path}")
                for row_number, row in enumerate(reader, 2):
                    accounting["read_level_rows_examined"] += 1
                    pair = (row.get("sample", ""), row.get("virus", ""))
                    if (
                        pair not in pairs
                        or row.get("mapping_mode") != "exact"
                        or row.get("virus_assignment") != "assigned"
                        or row.get("strand") != "antisense"
                    ):
                        continue
                    try:
                        length = int(row["length"])
                        abundance = float(row["count"])
                    except (TypeError, ValueError) as exc:
                        raise ValueError(f"invalid length/count in {path}:{row_number}") from exc
                    if length not in LENGTHS:
                        continue
                    sequence = row["sequence"].upper().replace("U", "T")
                    if len(sequence) != length or any(base not in DNA for base in sequence):
                        accounting["excluded_invalid_sequence_rows"] += 1
                        continue
                    if not math.isfinite(abundance) or abundance <= 0:
                        raise ValueError(f"non-positive/non-finite abundance in {path}:{row_number}")
                    output[(pair[0], pair[1], length)][sequence] += abundance
                    accounting["retained_observed_rows"] += 1
        except (gzip.BadGzipFile, EOFError, OSError) as exc:
            raise ValueError(f"corrupt frozen read-level gzip: {path}") from exc
    return output, dict(accounting)


def _format_number(value: object) -> object:
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def reconstruct_training_universe(
    legacy_core: Path,
    training_output: Path | None = None,
) -> dict[str, object]:
    pairs = load_primary_pairs(legacy_core)
    observed, scan_accounting = aggregate_observed_antisense(legacy_core, pairs)
    counts = Counter()
    represented_by_length = Counter()
    opportunities_by_length = Counter()
    supported_abundance = 0.0
    outside_species = 0
    outside_abundance = 0.0

    fieldnames = [
        "sample", "analysis_unit", "biological_virus", "candidate_length_nt",
        "sequence_dna", "represented", "abundance",
        *BASE_FEATURE_NAMES,
    ]
    handle = None
    writer = None
    if training_output is not None:
        training_output.parent.mkdir(parents=True, exist_ok=True)
        handle = gzip.open(training_output, "wt", newline="", encoding="utf-8")
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
    try:
        for (sample, unit), metadata in sorted(pairs.items()):
            fasta = legacy_core / "references/consensus" / f"{sample}.{unit}.final.background_masked.fa"
            if not fasta.is_file():
                raise FileNotFoundError(f"missing frozen background FASTA: {fasta}")
            records = parse_fasta(fasta)
            for length in LENGTHS:
                key = (sample, unit, length)
                opportunities = supported_antisense_sequences(records, length)
                abundance_by_sequence = observed.get(key, Counter())
                opportunities_by_length[length] += len(opportunities)
                matched = set(abundance_by_sequence).intersection(opportunities)
                represented_by_length[length] += len(matched)
                supported_abundance += sum(abundance_by_sequence[sequence] for sequence in matched)
                outside = set(abundance_by_sequence).difference(opportunities)
                outside_species += len(outside)
                outside_abundance += sum(abundance_by_sequence[sequence] for sequence in outside)
                if writer is not None:
                    for sequence in sorted(opportunities):
                        abundance = float(abundance_by_sequence.get(sequence, 0.0))
                        features = encode_predictors(guide_predictors(sequence))
                        writer.writerow({
                            "sample": sample,
                            "analysis_unit": unit,
                            "biological_virus": metadata["biological_virus"],
                            "candidate_length_nt": length,
                            "sequence_dna": sequence,
                            "represented": int(abundance > 0),
                            "abundance": _format_number(abundance),
                            **features,
                        })
    finally:
        if handle is not None:
            handle.close()

    counts.update({
        "primary_samples": len({sample for sample, _unit in pairs}),
        "sample_virus_units": len(pairs),
        "sample_virus_length_groups": len(pairs) * len(LENGTHS),
        "opportunities_23nt": opportunities_by_length[23],
        "opportunities_24nt": opportunities_by_length[24],
        "opportunities_total": sum(opportunities_by_length.values()),
        "represented_23nt": represented_by_length[23],
        "represented_24nt": represented_by_length[24],
        "represented_total": sum(represented_by_length.values()),
        "supported_abundance": supported_abundance,
        "outside_background_species": outside_species,
        "outside_background_abundance": outside_abundance,
    })
    counts.update(scan_accounting)
    return dict(counts)


def validate_frozen_accounting(accounting: Mapping[str, object]) -> None:
    failures = []
    for key, expected in EXPECTED_ACCOUNTING.items():
        observed = accounting.get(key)
        if observed is None or float(observed) != float(expected):
            failures.append(f"{key}: observed={observed!r}, expected={expected!r}")
    if failures:
        raise ValueError("Stage 09A frozen accounting regression failed: " + "; ".join(failures))


def sample_aware_weights(rows: Sequence[Mapping[str, object]]) -> list[float]:
    """Exact positive-sequence weights for one training fold."""
    if not rows:
        raise ValueError("cannot calculate weights for an empty training fold")
    group_keys = [
        (str(row["sample"]), str(row["analysis_unit"]), int(row["candidate_length_nt"]))
        for row in rows
    ]
    groups_by_sample: dict[str, set[tuple[str, str, int]]] = defaultdict(set)
    group_sizes = Counter(group_keys)
    for key in group_keys:
        groups_by_sample[key[0]].add(key)
    raw = [1.0 / (len(groups_by_sample[key[0]]) * group_sizes[key]) for key in group_keys]
    factor = len(rows) / sum(raw)
    return [value * factor for value in raw]


def fit_scaling(rows: Sequence[Mapping[str, object]], feature_names: Sequence[str] = BASE_FEATURE_NAMES) -> Scaling:
    if not rows:
        raise ValueError("cannot fit scaling on an empty training fold")
    means: dict[str, float] = {}
    sds: dict[str, float] = {}
    retained, omitted = [], []
    for name in feature_names:
        values = [float(row[name]) for row in rows]
        mean = statistics.fmean(values)
        sd = statistics.stdev(values) if len(values) > 1 else 0.0
        means[name], sds[name] = mean, sd
        (retained if sd > 0 else omitted).append(name)
    return Scaling(means, sds, tuple(retained), tuple(omitted))


def apply_scaling(row: Mapping[str, object], scaling: Scaling) -> dict[str, float]:
    return {name: (float(row[name]) - scaling.means[name]) / scaling.sds[name] for name in scaling.retained}


def weighted_within_group_center(
    values: Sequence[Sequence[float]],
    responses: Sequence[float],
    weights: Sequence[float],
    groups: Sequence[object],
) -> tuple[list[list[float]], list[float]]:
    if not (len(values) == len(responses) == len(weights) == len(groups)):
        raise ValueError("centering arrays must have identical lengths")
    if not values:
        raise ValueError("cannot center an empty training fold")
    width = len(values[0])
    members: dict[object, list[int]] = defaultdict(list)
    for index, group in enumerate(groups):
        members[group].append(index)
    centered_x = [[0.0] * width for _ in values]
    centered_y = [0.0] * len(values)
    for indexes in members.values():
        total_weight = sum(weights[index] for index in indexes)
        mean_y = sum(weights[index] * responses[index] for index in indexes) / total_weight
        mean_x = [
            sum(weights[index] * float(values[index][column]) for index in indexes) / total_weight
            for column in range(width)
        ]
        for index in indexes:
            centered_y[index] = responses[index] - mean_y
            centered_x[index] = [float(values[index][column]) - mean_x[column] for column in range(width)]
    return centered_x, centered_y


def structure_feature_names(structure: str, retained: Sequence[str], length: int | None = None) -> tuple[str, ...]:
    if structure == "A_shared":
        return tuple(retained)
    if structure == "B_shared_plus_length_interactions":
        return tuple(retained) + tuple(f"{name}_x_24nt" for name in retained)
    if structure == "C_separate_23_24":
        if length not in LENGTHS:
            raise ValueError("separate structure requires length 23 or 24")
        return tuple(retained)
    raise ValueError(f"unknown Stage 09A model structure: {structure}")


def structure_row(scaled: Mapping[str, float], structure: str, length: int) -> list[float]:
    names = tuple(scaled)
    values = [float(scaled[name]) for name in names]
    if structure == "A_shared" or structure == "C_separate_23_24":
        return values
    if structure == "B_shared_plus_length_interactions":
        return values + [value * int(length == 24) for value in values]
    raise ValueError(f"unknown Stage 09A model structure: {structure}")


def average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = ((start + 1) + end) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def favourable_percentiles(values: Sequence[float]) -> list[float]:
    if not values:
        return []
    ranks = average_ranks(values)
    return [(rank - 0.5) / len(values) for rank in ranks]


def top10_abundance_metrics(scores: Sequence[float], abundance: Sequence[float]) -> tuple[float | None, float | None]:
    if len(scores) != len(abundance) or not scores:
        raise ValueError("scores and abundance must have equal non-zero length")
    total = sum(abundance)
    if total <= 0:
        return None, None
    selected_n = math.ceil(0.10 * len(scores))
    order = sorted(range(len(scores)), key=lambda index: (-scores[index], index))[:selected_n]
    share = sum(abundance[index] for index in order) / total
    return share, share / (selected_n / len(scores))


def choose_model_configuration(rows: Sequence[Mapping[str, object]], tolerance: float = 1e-6) -> Mapping[str, object]:
    if not rows:
        raise ValueError("no model configurations to select")
    preference = {name: index for index, name in enumerate(STRUCTURES)}
    remaining = list(rows)
    best_rho = max(float(row["selection_score_rho"]) for row in remaining)
    remaining = [row for row in remaining if best_rho - float(row["selection_score_rho"]) <= tolerance]
    best_top10 = max(float(row["selection_score_top10"]) for row in remaining)
    remaining = [row for row in remaining if best_top10 - float(row["selection_score_top10"]) <= tolerance]
    return sorted(
        remaining,
        key=lambda row: (
            -float(row["alpha"]),
            -float(row["l1_ratio"]),
            preference[str(row["model_structure"])],
        ),
    )[0]


def assert_cv_partition(train_rows: Sequence[Mapping[str, object]], heldout_rows: Sequence[Mapping[str, object]], field: str) -> None:
    overlap = {str(row[field]) for row in train_rows}.intersection(str(row[field]) for row in heldout_rows)
    if overlap:
        raise ValueError(f"cross-validation leakage in {field}: {sorted(overlap)}")


def assert_no_stage08_feature_leakage(columns: Iterable[str]) -> None:
    leakage = sorted(set(columns).intersection(STAGE08_FORBIDDEN))
    if leakage:
        raise ValueError(f"Stage 08 features are forbidden in Stage 09A Layer 1: {leakage}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_rscript(environment: Mapping[str, str] | None = None) -> Path:
    """Resolve Rscript from the active reproducible Stage 09A environment."""
    environment = os.environ if environment is None else environment
    prefix = environment.get("CONDA_PREFIX", "")
    if not prefix:
        raise RuntimeError("Stage 09A requires an active conda environment with CONDA_PREFIX")
    rscript = Path(prefix) / "bin/Rscript"
    if not rscript.is_file():
        raise FileNotFoundError(f"Stage 09A environment lacks Rscript: {rscript}")
    return rscript


def prepare_candidates(path: Path, output: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(f"missing frozen Stage 06 candidate table: {path}")
    rows = _read_tsv(path)
    if not rows:
        raise ValueError("Stage 06 candidate table is empty")
    required = {
        "target_id", "candidate_id", "candidate_length_nt", "target_sequence_rna",
        "antisense_guide_sequence_rna",
    }
    if not required.issubset(rows[0]):
        raise ValueError("Stage 06 candidate table lacks required Stage 09A columns")
    assert_no_stage08_feature_leakage(rows[0])
    fieldnames = list(rows[0]) + list(BASE_FEATURE_NAMES)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        seen: set[str] = set()
        for row in rows:
            length = int(row["candidate_length_nt"])
            if length not in LENGTHS:
                raise ValueError(f"unsupported Stage 09A candidate length: {length}")
            guide = row["antisense_guide_sequence_rna"].upper()
            expected = reverse_complement_rna(row["target_sequence_rna"])
            if guide != expected or len(guide) != length:
                raise ValueError(f"candidate guide orientation mismatch: {row['candidate_id']}")
            if row["candidate_id"] in seen:
                raise ValueError(f"duplicate candidate_id: {row['candidate_id']}")
            seen.add(row["candidate_id"])
            writer.writerow({**row, **encode_predictors(guide_predictors(guide))})
    return len(rows)


def write_accounting(path: Path, accounting: Mapping[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in sorted(accounting):
            writer.writerow([key, _format_number(accounting[key])])


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--model-script", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    for path, label in ((args.solver, "approved glmnet adapter"), (args.model_script, "Stage 09A R model")):
        if not path.is_file():
            raise FileNotFoundError(f"missing {label}: {path}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="stage09a_", dir=args.output_root.parent) as temporary:
        temp = Path(temporary)
        training = temp / "training_universe.tsv.gz"
        candidate_features = temp / "candidate_features.tsv"
        accounting_path = temp / "accounting.tsv"
        accounting = reconstruct_training_universe(args.legacy_core, training)
        validate_frozen_accounting(accounting)
        candidate_count = prepare_candidates(args.candidates, candidate_features)
        accounting["candidate_rows"] = candidate_count
        accounting["stage06_candidates_sha256"] = sha256_file(args.candidates)
        write_accounting(accounting_path, accounting)
        command = [
            str(resolve_rscript()), str(args.model_script),
            "--training", str(training),
            "--candidates", str(candidate_features),
            "--accounting", str(accounting_path),
            "--solver", str(args.solver),
            "--output-root", str(args.output_root),
        ]
        subprocess.run(command, check=True)
    expected = [
        "layer1_model_coefficients.tsv", "layer1_model_preprocessing.tsv",
        "layer1_model_selection.tsv", "layer1_cv_by_group.tsv",
        "layer1_cv_summary_23nt.tsv", "layer1_cv_summary_24nt.tsv",
        "layer1_representation_diagnostic.tsv", "layer1_architecture_benchmarks.tsv",
        "candidate_layer1.tsv",
    ]
    missing = [name for name in expected if not (args.output_root / name).is_file()]
    if missing:
        raise RuntimeError(f"Stage 09A R model did not create required outputs: {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
