#!/usr/bin/env python3
"""Canonical Stage 06: generic exhaustive transcript-target enumeration."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from pathlib import Path
from typing import Iterable, Sequence


TARGET_MANIFEST_COLUMNS = [
    "target_id",
    "transcript_id",
    "display_name",
    "organism",
    "molecule_type",
    "fasta_path",
    "fasta_record_id",
    "annotation_path",
    "expected_length_nt",
    "sequence_sha256_uppercase_dna",
    "candidate_lengths_nt",
    "source_database",
    "source_accession_version",
]

ANNOTATION_COLUMNS = [
    "transcript_id",
    "region_label",
    "start_1based",
    "end_1based",
    "coordinate_system",
]

CANDIDATE_COLUMNS = [
    "target_id",
    "transcript_id",
    "display_name",
    "organism",
    "candidate_id",
    "candidate_length_nt",
    "start_1based",
    "end_1based",
    "target_sequence_dna",
    "target_sequence_rna",
    "antisense_guide_sequence_rna",
    "annotation_status",
    "start_region",
    "end_region",
    "overlap_regions",
    "crosses_annotation_boundary",
]

REFERENCE_SUMMARY_COLUMNS = [
    "target_id",
    "transcript_id",
    "display_name",
    "organism",
    "molecule_type",
    "fasta_path",
    "fasta_record_id",
    "annotation_path",
    "normalized_length_nt",
    "sequence_sha256_uppercase_dna",
    "candidate_lengths_nt",
    "annotation_status",
    "candidate_count",
]

PROVENANCE_COLUMNS = TARGET_MANIFEST_COLUMNS + [
    "observed_length_nt",
    "observed_sequence_sha256_uppercase_dna",
    "annotation_status",
    "candidate_count",
    "validation_status",
]

QC_COLUMNS = ["status", "check", "target_id", "candidate_length_nt", "value", "details"]
COORDINATE_SYSTEM = "1-based inclusive transcript coordinates"
NA_VALUES = {"", "NA", "N/A", "NONE", "."}


class Stage06Error(ValueError):
    """A deterministic Stage 06 validation error."""


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames or []
            return header, list(reader)
    except OSError as exc:
        raise Stage06Error(f"Cannot read TSV {path}: {exc}") from exc


def _resolve_path(manifest_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    repository_relative = Path.cwd() / path
    manifest_relative = manifest_path.parent / path
    if repository_relative.exists():
        return repository_relative.resolve()
    return manifest_relative.resolve()


def load_target_manifest(path: str | Path) -> list[dict[str, str]]:
    """Load and validate the generic target registry schema and identifiers."""
    manifest_path = Path(path).resolve()
    header, rows = _read_tsv(manifest_path)
    if header != TARGET_MANIFEST_COLUMNS:
        raise Stage06Error(
            "target manifest schema mismatch: expected "
            f"{TARGET_MANIFEST_COLUMNS}, observed {header}"
        )
    if not rows:
        raise Stage06Error("target manifest contains no targets")

    seen_target_ids: set[str] = set()
    loaded: list[dict[str, str]] = []
    for line_number, source_row in enumerate(rows, start=2):
        row = {key: (value or "").strip() for key, value in source_row.items()}
        missing = [key for key in TARGET_MANIFEST_COLUMNS if not row[key]]
        if missing:
            raise Stage06Error(f"manifest line {line_number} has empty fields: {missing}")
        if row["target_id"] in seen_target_ids:
            raise Stage06Error(f"duplicate target_id: {row['target_id']}")
        seen_target_ids.add(row["target_id"])
        row["_manifest_path"] = str(manifest_path)
        row["_fasta_resolved"] = str(_resolve_path(manifest_path, row["fasta_path"]))
        annotation = row["annotation_path"].upper()
        row["_annotation_resolved"] = (
            "" if annotation in NA_VALUES else str(_resolve_path(manifest_path, row["annotation_path"]))
        )
        loaded.append(row)
    return loaded


def normalize_sequence(sequence: str) -> str:
    """Normalize DNA/RNA sequence text to uppercase unambiguous DNA."""
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    if not normalized:
        raise Stage06Error("transcript sequence is empty")
    unexpected = sorted(set(normalized) - set("ACGT"))
    if unexpected:
        raise Stage06Error(f"transcript contains ambiguous/unexpected bases: {unexpected}")
    return normalized


def load_transcript_sequence(fasta_path: str | Path, fasta_record_id: str) -> str:
    """Select one exact FASTA record and return normalized uppercase DNA."""
    path = Path(fasta_path)
    records: dict[str, list[str]] = {}
    current_id: str | None = None
    try:
        with path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    record_id = line[1:].split(maxsplit=1)[0]
                    if not record_id:
                        raise Stage06Error(f"empty FASTA identifier at {path}:{line_number}")
                    if record_id in records:
                        raise Stage06Error(f"duplicate FASTA record identifier: {record_id}")
                    records[record_id] = []
                    current_id = record_id
                elif current_id is None:
                    raise Stage06Error(f"sequence occurs before a FASTA header at {path}:{line_number}")
                else:
                    records[current_id].append(line)
    except OSError as exc:
        raise Stage06Error(f"Cannot read FASTA {path}: {exc}") from exc

    if fasta_record_id not in records:
        raise Stage06Error(
            f"requested FASTA record {fasta_record_id!r} not found in {path}; "
            f"available records: {sorted(records)}"
        )
    return normalize_sequence("".join(records[fasta_record_id]))


def parse_candidate_lengths(value: str, transcript_length: int) -> list[int]:
    """Parse a comma-separated list of distinct positive candidate lengths."""
    parsed: list[int] = []
    for token in value.split(","):
        token = token.strip()
        try:
            length = int(token)
        except ValueError as exc:
            raise Stage06Error(f"invalid candidate length: {token!r}") from exc
        if length <= 0:
            raise Stage06Error(f"candidate length must be positive: {length}")
        if length > transcript_length:
            raise Stage06Error(
                f"candidate length {length} exceeds transcript length {transcript_length}"
            )
        if length in parsed:
            raise Stage06Error(f"duplicate candidate length: {length}")
        parsed.append(length)
    if not parsed:
        raise Stage06Error("no candidate lengths requested")
    return parsed


def sequence_sha256(sequence_dna: str) -> str:
    return hashlib.sha256(sequence_dna.encode("ascii")).hexdigest()


def load_transcript_regions(
    annotation_path: str | Path | None,
    transcript_id: str,
    transcript_length: int,
) -> tuple[list[dict[str, int | str]] | None, str]:
    """Load optional, non-overlapping transcript-coordinate annotations."""
    if annotation_path is None or str(annotation_path).strip().upper() in NA_VALUES:
        return None, "unavailable"

    path = Path(annotation_path)
    header, all_rows = _read_tsv(path)
    if header != ANNOTATION_COLUMNS:
        raise Stage06Error(
            f"annotation schema mismatch for {path}: expected {ANNOTATION_COLUMNS}, observed {header}"
        )
    selected = [row for row in all_rows if (row["transcript_id"] or "").strip() == transcript_id]
    if not selected:
        raise Stage06Error(f"annotation has no rows for transcript_id {transcript_id!r}")

    regions: list[dict[str, int | str]] = []
    for row in selected:
        label = (row["region_label"] or "").strip()
        if not label:
            raise Stage06Error("annotation region_label must not be empty")
        if (row["coordinate_system"] or "").strip() != COORDINATE_SYSTEM:
            raise Stage06Error(
                f"unsupported annotation coordinate system: {row['coordinate_system']!r}"
            )
        try:
            start = int(row["start_1based"])
            end = int(row["end_1based"])
        except (TypeError, ValueError) as exc:
            raise Stage06Error(f"non-integer annotation coordinates: {row}") from exc
        if start < 1 or end < start or end > transcript_length:
            raise Stage06Error(
                f"annotation interval {label} [{start}, {end}] lies outside transcript "
                f"1..{transcript_length}"
            )
        regions.append({"region_label": label, "start_1based": start, "end_1based": end})

    regions.sort(key=lambda item: (int(item["start_1based"]), int(item["end_1based"])))
    for previous, current in zip(regions, regions[1:]):
        if int(current["start_1based"]) <= int(previous["end_1based"]):
            raise Stage06Error(
                "annotation intervals overlap: "
                f"{previous['region_label']} and {current['region_label']}"
            )

    complete = (
        int(regions[0]["start_1based"]) == 1
        and int(regions[-1]["end_1based"]) == transcript_length
        and all(
            int(right["start_1based"]) == int(left["end_1based"]) + 1
            for left, right in zip(regions, regions[1:])
        )
    )
    return regions, "complete" if complete else "partial"


def validate_target(
    target: dict[str, str],
    sequence_dna: str,
) -> list[int]:
    """Validate registry identity metadata and requested lengths."""
    try:
        expected_length = int(target["expected_length_nt"])
    except ValueError as exc:
        raise Stage06Error(
            f"invalid expected_length_nt for {target['target_id']}: {target['expected_length_nt']!r}"
        ) from exc
    if expected_length <= 0:
        raise Stage06Error(f"expected_length_nt must be positive for {target['target_id']}")
    if len(sequence_dna) != expected_length:
        raise Stage06Error(
            f"length mismatch for {target['target_id']}: expected {expected_length}, "
            f"observed {len(sequence_dna)}"
        )
    observed_hash = sequence_sha256(sequence_dna)
    expected_hash = target["sequence_sha256_uppercase_dna"].lower()
    if observed_hash != expected_hash:
        raise Stage06Error(
            f"SHA-256 mismatch for {target['target_id']}: expected {expected_hash}, "
            f"observed {observed_hash}"
        )
    return parse_candidate_lengths(target["candidate_lengths_nt"], len(sequence_dna))


def reverse_complement_rna(sequence_rna: str) -> str:
    """Return an RNA reverse complement in biological 5' to 3' orientation."""
    sequence = "".join(sequence_rna.split()).upper().replace("T", "U")
    unexpected = sorted(set(sequence) - set("ACGU"))
    if unexpected:
        raise Stage06Error(f"RNA contains ambiguous/unexpected bases: {unexpected}")
    return sequence.translate(str.maketrans("ACGU", "UGCA"))[::-1]


def _region_at(position: int, regions: Sequence[dict[str, int | str]]) -> str:
    for region in regions:
        if int(region["start_1based"]) <= position <= int(region["end_1based"]):
            return str(region["region_label"])
    return "unannotated"


def annotate_candidate(
    start_1based: int,
    end_1based: int,
    regions: Sequence[dict[str, int | str]] | None,
    annotation_status: str,
) -> dict[str, str]:
    """Attach optional region labels without filtering boundary-spanning candidates."""
    if regions is None:
        return {
            "annotation_status": "unavailable",
            "start_region": "NA",
            "end_region": "NA",
            "overlap_regions": "NA",
            "crosses_annotation_boundary": "NA",
        }

    labels: list[str] = []
    for position in range(start_1based, end_1based + 1):
        label = _region_at(position, regions)
        if not labels or labels[-1] != label:
            labels.append(label)
    return {
        "annotation_status": annotation_status,
        "start_region": labels[0],
        "end_region": labels[-1],
        "overlap_regions": ";".join(labels),
        "crosses_annotation_boundary": "TRUE" if len(labels) > 1 else "FALSE",
    }


def enumerate_candidates(
    target: dict[str, str],
    sequence_dna: str,
    candidate_lengths: Iterable[int],
    regions: Sequence[dict[str, int | str]] | None,
    annotation_status: str,
) -> list[dict[str, str | int]]:
    """Enumerate every complete transcript interval for each requested length."""
    transcript_length = len(sequence_dna)
    coordinate_width = max(4, len(str(transcript_length)))
    candidates: list[dict[str, str | int]] = []
    for candidate_length in candidate_lengths:
        for start in range(1, transcript_length - candidate_length + 2):
            end = start + candidate_length - 1
            target_dna = sequence_dna[start - 1 : end]
            target_rna = target_dna.replace("T", "U")
            candidate_id = (
                f"{target['target_id']}__{candidate_length}nt__"
                f"{start:0{coordinate_width}d}_{end:0{coordinate_width}d}"
            )
            row: dict[str, str | int] = {
                "target_id": target["target_id"],
                "transcript_id": target["transcript_id"],
                "display_name": target["display_name"],
                "organism": target["organism"],
                "candidate_id": candidate_id,
                "candidate_length_nt": candidate_length,
                "start_1based": start,
                "end_1based": end,
                "target_sequence_dna": target_dna,
                "target_sequence_rna": target_rna,
                "antisense_guide_sequence_rna": reverse_complement_rna(target_rna),
            }
            row.update(annotate_candidate(start, end, regions, annotation_status))
            candidates.append(row)
    return candidates


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _qc_row(
    status: str,
    check: str,
    value: object,
    details: str,
    target_id: str = "ALL",
    candidate_length_nt: object = "NA",
) -> dict[str, object]:
    return {
        "status": status,
        "check": check,
        "target_id": target_id,
        "candidate_length_nt": candidate_length_nt,
        "value": value,
        "details": details,
    }


def run_stage06(target_manifest: str | Path, output_root: str | Path) -> dict[str, int]:
    """Validate all registered targets, enumerate candidates, and write Stage 06 outputs."""
    output_root = Path(output_root)
    targets = load_target_manifest(target_manifest)
    all_candidates: list[dict[str, str | int]] = []
    summaries: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    qc_rows = [
        _qc_row("PASS", "target_manifest_schema", len(targets), "required schema is exact"),
        _qc_row("PASS", "registered_targets", len(targets), "target IDs are non-empty and unique"),
    ]

    for target in targets:
        sequence = load_transcript_sequence(target["_fasta_resolved"], target["fasta_record_id"])
        lengths = validate_target(target, sequence)
        annotation_path = target["_annotation_resolved"] or None
        regions, annotation_status = load_transcript_regions(
            annotation_path, target["transcript_id"], len(sequence)
        )
        candidates = enumerate_candidates(
            target, sequence, lengths, regions, annotation_status
        )
        all_candidates.extend(candidates)
        observed_hash = sequence_sha256(sequence)
        expected_total = sum(len(sequence) - length + 1 for length in lengths)

        summary = {
            **target,
            "normalized_length_nt": len(sequence),
            "sequence_sha256_uppercase_dna": observed_hash,
            "annotation_status": annotation_status,
            "candidate_count": len(candidates),
        }
        summaries.append(summary)
        provenance.append(
            {
                **target,
                "observed_length_nt": len(sequence),
                "observed_sequence_sha256_uppercase_dna": observed_hash,
                "annotation_status": annotation_status,
                "candidate_count": len(candidates),
                "validation_status": "PASS",
            }
        )
        qc_rows.append(
            _qc_row(
                "PASS",
                "target_reference_identity",
                len(sequence),
                f"length and normalized sequence SHA-256 match registry; annotation={annotation_status}",
                target["target_id"],
            )
        )
        qc_rows.append(
            _qc_row(
                "PASS",
                "target_candidate_total",
                len(candidates),
                f"expected exhaustive count={expected_total}",
                target["target_id"],
            )
        )
        if annotation_status == "unavailable":
            qc_rows.append(
                _qc_row(
                    "INFO",
                    "annotation_availability",
                    annotation_status,
                    "annotation is optional; candidates retain NA annotation fields",
                    target["target_id"],
                )
            )
        else:
            qc_rows.append(
                _qc_row(
                    "PASS",
                    "annotation_validation",
                    annotation_status,
                    "coordinates are in bounds and non-overlapping; gaps are permitted",
                    target["target_id"],
                )
            )

        for length in lengths:
            stratum = [row for row in candidates if row["candidate_length_nt"] == length]
            expected = len(sequence) - length + 1
            interval_ok = (
                len(stratum) == expected
                and int(stratum[0]["start_1based"]) == 1
                and int(stratum[0]["end_1based"]) == length
                and int(stratum[-1]["start_1based"]) == expected
                and int(stratum[-1]["end_1based"]) == len(sequence)
                and all(
                    row["target_sequence_dna"]
                    == sequence[int(row["start_1based"]) - 1 : int(row["end_1based"])]
                    and row["antisense_guide_sequence_rna"]
                    == reverse_complement_rna(str(row["target_sequence_rna"]))
                    for row in stratum
                )
            )
            if not interval_ok:
                raise Stage06Error(
                    f"internal exhaustive-enumeration invariant failed for {target['target_id']} "
                    f"length {length}"
                )
            qc_rows.append(
                _qc_row(
                    "PASS",
                    "candidate_stratum_enumeration",
                    len(stratum),
                    "all legal starts retained; slicing and guide orientation validated",
                    target["target_id"],
                    length,
                )
            )

    candidate_ids = [str(row["candidate_id"]) for row in all_candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise Stage06Error("candidate IDs are not globally unique")
    qc_rows.extend(
        [
            _qc_row(
                "PASS",
                "candidate_id_global_uniqueness",
                len(candidate_ids),
                "all candidate IDs are globally unique",
            ),
            _qc_row(
                "PASS",
                "canonical_candidate_schema",
                len(CANDIDATE_COLUMNS),
                "candidate table uses the exact canonical Stage 06 columns",
            ),
            _qc_row(
                "PASS",
                "score_rank_absence",
                0,
                "no score or rank columns are present",
            ),
            _qc_row(
                "PASS",
                "total_candidates",
                len(all_candidates),
                "unfiltered exhaustive candidate count",
            ),
        ]
    )

    _write_tsv(output_root / "target_reference_summary.tsv", REFERENCE_SUMMARY_COLUMNS, summaries)
    _write_tsv(output_root / "target_candidates.tsv", CANDIDATE_COLUMNS, all_candidates)
    _write_tsv(output_root / "provenance/stage06_manifest.tsv", PROVENANCE_COLUMNS, provenance)
    _write_tsv(output_root / "qc/stage06_accounting.tsv", QC_COLUMNS, qc_rows)
    return {
        "targets": len(targets),
        "candidates": len(all_candidates),
        "pass": sum(row["status"] == "PASS" for row in qc_rows),
        "warn": sum(row["status"] == "WARN" for row in qc_rows),
        "fail": sum(row["status"] == "FAIL" for row in qc_rows),
        "info": sum(row["status"] == "INFO" for row in qc_rows),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        run_stage06(args.target_manifest, args.output_root)
    except Stage06Error as exc:
        _write_tsv(
            args.output_root / "qc/stage06_accounting.tsv",
            QC_COLUMNS,
            [_qc_row("FAIL", "stage06_execution", 1, str(exc))],
        )
        print(f"Stage 06 validation failed: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
