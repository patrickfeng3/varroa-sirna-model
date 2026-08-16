#!/usr/bin/env python3
"""Export canonical Stage 10A evidence to the compact Stage 11 web schema."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence


SCHEMA_VERSION = "stage11-web-v1"
PIPELINE_SPEC_VERSION = "0.22"
METRIC_DICTIONARY_VERSION = "0.22"
METRIC_MAP = {
    "layer1": "stage10_layer1_percentile",
    "layer2": "stage10_layer2_percentile",
    "layer3": "stage10_layer3_percentile",
    "total": "stage10_equal_layer_score",
}
EXPORTED_CANDIDATE_FIELDS = {
    "candidate_id",
    "candidate_length_nt",
    "target_start_1based",
    "target_end_1based",
    "layer1",
    "layer2",
    "layer3",
    "total",
}
FEATURE_LABELS = {
    "5_prime_UTR": "5UTR",
    "5UTR": "5UTR",
    "5'UTR": "5UTR",
    "CDS": "CDS",
    "3_prime_UTR": "3UTR",
    "3UTR": "3UTR",
    "3'UTR": "3UTR",
}


def _parse_tsv_text(text: str, label: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"Missing TSV header: {label}")
    rows = list(reader)
    if not rows:
        raise ValueError(f"Empty TSV: {label}")
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required Stage 11 input is missing: {path}")
    return _parse_tsv_text(path.read_text(encoding="utf-8"), str(path))


def read_stage10_once(path: Path) -> tuple[list[dict[str, str]], str]:
    if not path.is_file():
        raise FileNotFoundError(f"Required frozen Stage 10A input is missing: {path}")
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    return _parse_tsv_text(raw.decode("utf-8"), str(path)), digest


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_columns(rows: Sequence[dict], columns: Sequence[str], label: str) -> None:
    missing = sorted(set(columns) - set(rows[0]))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def normalize_dna(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    if not normalized or set(normalized) - set("ACGT"):
        raise ValueError("Transcript sequence must contain only A/C/G/T or U")
    return normalized


def load_fasta_record(path: Path, record_id: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"Target FASTA is missing: {path}")
    records: dict[str, list[str]] = {}
    current: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            if current in records:
                raise ValueError(f"Duplicate FASTA record ID: {current}")
            records[current] = []
        elif current is None:
            raise ValueError(f"FASTA sequence before header: {path}")
        else:
            records[current].append(line)
    if record_id not in records:
        raise ValueError(f"FASTA record {record_id!r} not found in {path}")
    return normalize_dna("".join(records[record_id]))


def load_target_metadata(manifest_path: Path, target_id: str, repo_root: Path) -> dict:
    rows = read_tsv(manifest_path)
    required = [
        "target_id",
        "transcript_id",
        "display_name",
        "fasta_path",
        "fasta_record_id",
        "annotation_path",
        "expected_length_nt",
        "sequence_sha256_uppercase_dna",
    ]
    require_columns(rows, required, "target manifest")
    matches = [row for row in rows if row["target_id"] == target_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one target-manifest row for {target_id}; found {len(matches)}")
    row = matches[0]
    fasta_path = repo_root / row["fasta_path"]
    sequence = load_fasta_record(fasta_path, row["fasta_record_id"])
    sequence_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    expected_length = int(row["expected_length_nt"])
    if len(sequence) != expected_length:
        raise ValueError(f"Transcript length mismatch for {target_id}: {len(sequence)} != {expected_length}")
    if sequence_hash != row["sequence_sha256_uppercase_dna"]:
        raise ValueError(f"Transcript SHA-256 mismatch for {target_id}")

    annotation_path = repo_root / row["annotation_path"]
    annotations_raw = read_tsv(annotation_path)
    require_columns(annotations_raw, ["transcript_id", "region_label", "start_1based", "end_1based"], "annotation")
    annotations = []
    for annotation in annotations_raw:
        if annotation["transcript_id"] != row["transcript_id"]:
            continue
        label = annotation["region_label"]
        if label not in FEATURE_LABELS:
            raise ValueError(f"Unsupported Stage 11 display feature label: {label}")
        start = int(annotation["start_1based"])
        end = int(annotation["end_1based"])
        if not 1 <= start <= end <= len(sequence):
            raise ValueError(f"Annotation outside transcript: {label} {start}-{end}")
        annotations.append({"feature": FEATURE_LABELS[label], "start_1based": start, "end_1based": end})
    annotations.sort(key=lambda item: (item["start_1based"], item["end_1based"], item["feature"]))
    for left, right in zip(annotations, annotations[1:]):
        if left["end_1based"] >= right["start_1based"]:
            raise ValueError("Stage 11 annotations overlap")
    if {item["feature"] for item in annotations} != {"5UTR", "CDS", "3UTR"}:
        raise ValueError("Stage 11 requires canonical 5UTR, CDS and 3UTR annotations")
    return {
        "target_id": target_id,
        "transcript_id": row["transcript_id"],
        "display_name": row["display_name"],
        "transcript_length_nt": len(sequence),
        "transcript_sequence": sequence,
        "transcript_sequence_sha256": sequence_hash,
        "annotations": annotations,
    }


def _git_commit(repo_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return result.stdout.strip() or "unavailable"
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def build_payload(
    stage10_path: Path,
    manifest_path: Path,
    target_id: str,
    repo_root: Path,
    *,
    source_git_commit: str | None = None,
    enforce_current_fixture: bool = False,
) -> tuple[dict, list[dict]]:
    stage10_rows, stage10_hash = read_stage10_once(stage10_path)
    required = [
        "candidate_id",
        "target_id",
        "candidate_length_nt",
        "start_1based",
        "end_1based",
        "target_sequence_dna",
        *METRIC_MAP.values(),
    ]
    require_columns(stage10_rows, required, "Stage 10A candidates")
    target = load_target_metadata(manifest_path, target_id, repo_root)
    selected = [row for row in stage10_rows if row["target_id"] == target_id]
    if not selected:
        raise ValueError(f"No Stage 10A candidates for target {target_id}")
    candidate_ids = [row["candidate_id"] for row in selected]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("Duplicate Stage 10A candidate IDs")

    candidates = []
    by_length: dict[int, list[int]] = defaultdict(list)
    for row in selected:
        candidate_id = row["candidate_id"]
        length = int(row["candidate_length_nt"])
        start = int(row["start_1based"])
        end = int(row["end_1based"])
        if end != start + length - 1 or not 1 <= start <= end <= target["transcript_length_nt"]:
            raise ValueError(f"Invalid Stage 10A candidate coordinates: {candidate_id}")
        if row["target_sequence_dna"] != target["transcript_sequence"][start - 1 : end]:
            raise ValueError(f"Stage 10A target slice mismatch: {candidate_id}")
        scores = {}
        for web_name, source_name in METRIC_MAP.items():
            value = float(row[source_name])
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"Invalid Stage 10A score {source_name} for {candidate_id}")
            scores[web_name] = value
        candidates.append(
            {
                "candidate_id": candidate_id,
                "candidate_length_nt": length,
                "target_start_1based": start,
                "target_end_1based": end,
                **scores,
            }
        )
        by_length[length].append(start)

    supported_lengths = sorted(by_length)
    for length, starts in by_length.items():
        expected = list(range(1, target["transcript_length_nt"] - length + 2))
        if sorted(starts) != expected:
            raise ValueError(f"Stage 10A candidates are not exhaustive for guide length {length}")
    counts = Counter(candidate["candidate_length_nt"] for candidate in candidates)
    if enforce_current_fixture and (len(candidates) != 1375 or counts != {23: 688, 24: 687}):
        raise ValueError(f"Current Stage 11 fixture failed: total={len(candidates)}, lengths={dict(counts)}")

    candidates.sort(key=lambda item: (item["candidate_length_nt"], item["target_start_1based"], item["candidate_id"]))
    payload = {
        "schema_version": SCHEMA_VERSION,
        "target": target,
        "supported_guide_lengths": supported_lengths,
        "metrics": dict(METRIC_MAP),
        "candidates": candidates,
        "provenance": {
            "pipeline_spec_version": PIPELINE_SPEC_VERSION,
            "metric_dictionary_version": METRIC_DICTIONARY_VERSION,
            "source_stage10_git_commit": source_git_commit if source_git_commit is not None else _git_commit(repo_root),
            "source_stage10_file_sha256": stage10_hash,
        },
    }

    exported = {candidate["candidate_id"]: candidate for candidate in candidates}
    exact_scores = True
    for row in selected:
        candidate = exported[row["candidate_id"]]
        exact_scores = exact_scores and all(candidate[name] == float(row[source]) for name, source in METRIC_MAP.items())
    checks = [
        ("schema_version", payload["schema_version"] == SCHEMA_VERSION, payload["schema_version"], "canonical web schema"),
        ("candidate_ids_preserved", set(exported) == set(candidate_ids), len(exported), "exact one-to-one Stage 10A export"),
        ("candidate_scores_preserved", exact_scores, int(exact_scores), "L1/L2/L3/total exact"),
        ("candidate_count_total", len(candidates) == 1375, len(candidates), "current Vd-CHIBIN fixture"),
        ("candidate_count_23nt", counts.get(23) == 688, counts.get(23, 0), "current Vd-CHIBIN fixture"),
        ("candidate_count_24nt", counts.get(24) == 687, counts.get(24, 0), "current Vd-CHIBIN fixture"),
        ("transcript_length", target["transcript_length_nt"] == 710, target["transcript_length_nt"], "current Vd-CHIBIN fixture"),
        ("transcript_sha256", target["transcript_sequence_sha256"] == "4a0d25aa05b269a118ed1b952dca63ccd1c0a7978fc42295faf3bf650e43ea42", target["transcript_sequence_sha256"], "canonical target manifest"),
        ("annotations_present", {item["feature"] for item in target["annotations"]} == {"5UTR", "CDS", "3UTR"}, len(target["annotations"]), "canonical transcript metadata"),
        ("minimal_candidate_schema", all(set(item) == EXPORTED_CANDIDATE_FIELDS for item in candidates), len(EXPORTED_CANDIDATE_FIELDS), "no Pareto/minimum/raw upstream fields"),
        ("supported_guide_lengths", supported_lengths == [23, 24], ",".join(map(str, supported_lengths)), "never mixed in region score"),
    ]
    qc = [
        {"check": name, "status": "PASS" if passed else "FAIL", "value": value, "detail": detail}
        for name, passed, value, detail in checks
    ]
    return payload, qc


def feature_for_position(annotations: Sequence[dict], position: int) -> str | None:
    for annotation in annotations:
        if annotation["start_1based"] <= position <= annotation["end_1based"]:
            return annotation["feature"]
    return None


def compute_regions(payload: Mapping, guide_length: int, region_length: int, metric: str) -> list[dict]:
    target = payload["target"]
    transcript_length = int(target["transcript_length_nt"])
    if guide_length not in payload["supported_guide_lengths"] or metric not in METRIC_MAP:
        return []
    if not isinstance(region_length, int) or region_length < guide_length or region_length > transcript_length:
        return []
    selected = sorted(
        (candidate for candidate in payload["candidates"] if candidate["candidate_length_nt"] == guide_length),
        key=lambda item: item["target_start_1based"],
    )
    expected_starts = list(range(1, transcript_length - guide_length + 2))
    if [candidate["target_start_1based"] for candidate in selected] != expected_starts:
        raise ValueError(f"Incomplete candidate starts for guide length {guide_length}")
    values = [float(candidate[metric]) for candidate in selected]
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    contained_count = region_length - guide_length + 1
    regions = []
    sequence = target["transcript_sequence"]
    for start in range(1, transcript_length - region_length + 2):
        end = start + region_length - 1
        first_index = start - 1
        mean_score = (
            values[first_index]
            if contained_count == 1
            else (prefix[first_index + contained_count] - prefix[first_index]) / contained_count
        )
        regions.append(
            {
                "stage11_region_start_1based": start,
                "stage11_region_end_1based": end,
                "stage11_region_start_feature": feature_for_position(target["annotations"], start),
                "stage11_guide_length_nt": guide_length,
                "stage11_region_length_nt": region_length,
                "stage11_metric_mode": metric,
                "stage11_region_mean_score": mean_score,
                "stage11_region_n_contained_windows": contained_count,
                "stage11_region_sequence": sequence[start - 1 : end],
            }
        )
    return regions


def top_bottom_by_feature(
    regions: Sequence[dict], features: Sequence[str] = ("5UTR", "CDS", "3UTR"), limit: int = 5
) -> dict[str, dict[str, list[dict | None]]]:
    output = {}
    for feature in features:
        matching = [region for region in regions if region["stage11_region_start_feature"] == feature]
        top = sorted(matching, key=lambda item: (-item["stage11_region_mean_score"], item["stage11_region_start_1based"]))[:limit]
        bottom = sorted(matching, key=lambda item: (item["stage11_region_mean_score"], item["stage11_region_start_1based"]))[:limit]
        output[feature] = {
            "top": top + [None] * (limit - len(top)),
            "bottom": bottom + [None] * (limit - len(bottom)),
        }
    return output


def validate_url_state(
    query: Mapping[str, str], supported_lengths: Sequence[int], transcript_length: int
) -> dict[str, int | str]:
    default_guide = min(supported_lengths)
    try:
        supplied_guide = int(query.get("guide", ""))
    except (TypeError, ValueError):
        supplied_guide = default_guide
    guide = supplied_guide if supplied_guide in supported_lengths else default_guide
    default_region = min(max(96, guide), transcript_length)
    try:
        supplied_region = int(query.get("region", ""))
    except (TypeError, ValueError):
        supplied_region = default_region
    region = supplied_region if guide <= supplied_region <= transcript_length else default_region
    supplied_metric = query.get("metric", "total")
    metric = supplied_metric if supplied_metric in METRIC_MAP else "total"
    return {"guide": guide, "region": region, "metric": metric}


def compact_json(payload: Mapping) -> str:
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"


def compact_javascript(payload: Mapping) -> str:
    """Equivalent static payload for browsers opened directly with file://."""
    return "window.STAGE11_PAYLOAD=" + compact_json(payload).rstrip("\n") + ";\n"


def run(args: argparse.Namespace) -> None:
    repo_root = args.repo_root.resolve()
    payload, qc = build_payload(
        args.stage10,
        args.target_manifest,
        args.target_id,
        repo_root,
        enforce_current_fixture=True,
    )
    failures = [row for row in qc if row["status"] == "FAIL"]
    if failures:
        raise ValueError("Stage 11 export QC failed: " + "; ".join(row["check"] for row in failures))
    args.web_data.parent.mkdir(parents=True, exist_ok=True)
    args.web_data.write_text(compact_json(payload), encoding="utf-8")
    args.web_data_js.parent.mkdir(parents=True, exist_ok=True)
    args.web_data_js.write_text(compact_javascript(payload), encoding="utf-8")
    write_tsv(args.qc, qc, ["check", "status", "value", "detail"])
    parameters = [
        {"parameter": "schema_version", "value": SCHEMA_VERSION, "detail": "static web-data contract"},
        {"parameter": "metric_mapping", "value": "layer1:L1;layer2:L2;layer3:L3;total:equal_layer", "detail": "Stage 10A values only"},
        {"parameter": "region_score", "value": "unweighted_mean_fully_contained_selected_length_windows", "detail": "no Pareto/minimum penalty"},
        {"parameter": "feature_assignment", "value": "region_start_coordinate", "detail": "boundary-crossing regions retained"},
        {"parameter": "smoothing", "value": "none_beyond_requested_region_mean", "detail": "R-L+1 windows"},
        {"parameter": "stage10_source", "value": str(args.stage10), "detail": "read-only; not a DAG dependency"},
        {"parameter": "target_manifest", "value": str(args.target_manifest), "detail": "canonical transcript metadata"},
    ]
    write_tsv(args.parameters, parameters, ["parameter", "value", "detail"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage10", type=Path, required=True)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--target-id", required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--web-data", type=Path, required=True)
    parser.add_argument("--web-data-js", type=Path, required=True)
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
