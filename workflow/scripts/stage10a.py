#!/usr/bin/env python3
"""Canonical Stage 10A individual-window evidence integration."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence


IDENTITY_COLUMNS = [
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

CANDIDATE_COLUMNS = IDENTITY_COLUMNS + [
    "layer1_accumulation_percentile",
    "layer2_reference_score",
    "layer3_reference_score",
    "stage10_layer1_percentile",
    "stage10_layer2_percentile",
    "stage10_layer3_percentile",
    "stage10_equal_layer_score",
    "stage10_equal_layer_rank",
    "stage10_equal_layer_percentile",
    "stage10_pareto_front",
    "stage10_minimum_layer_score",
]

FORBIDDEN_STAGE08_RAW_COLUMNS = {
    "target_whole_p_unpaired",
    "target_whole_p_unpaired_w100_l80",
    "target_whole_p_unpaired_w200_l150",
    "target_seed_g2_8_p_unpaired",
    "target_seed_g2_8_p_unpaired_w100_l80",
    "target_seed_g2_8_p_unpaired_w200_l150",
    "guide_5p_terminal_dg_4bp",
    "passenger_5p_terminal_dg_4bp",
    "asymmetry_ddg_4bp",
    "guide_5p_terminal_dg_5bp",
    "passenger_5p_terminal_dg_5bp",
    "asymmetry_ddg_5bp",
    "guide_self_fold_mfe_kcal_mol",
    "guide_self_fold_structure",
}

FORBIDDEN_OUTPUT_COLUMNS = {
    "efficacy_probability",
    "overall_score",
    "overall_rank",
    "region_score",
    "selected_region",
    "construct_score",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required frozen Stage 09 input is missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Stage 09 input is empty: {path}")
    return rows


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_columns(rows: Sequence[dict], columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(rows[0]))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def average_ranks(values: Sequence[float], *, descending: bool = False) -> list[float]:
    """Return one-based average ranks, preserving exact ties."""
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda index: values[index], reverse=descending)
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        value = values[order[position]]
        while end < len(order) and values[order[end]] == value:
            end += 1
        average_rank = ((position + 1) + end) / 2.0
        for index in order[position:end]:
            ranks[index] = average_rank
        position = end
    return ranks


def favourable_percentiles(values: Sequence[float]) -> list[float]:
    if not values:
        return []
    n = len(values)
    return [(rank - 0.5) / n for rank in average_ranks(values)]


def pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    if len(x) != len(y) or len(x) < 2:
        return None
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    centered_x = [value - mean_x for value in x]
    centered_y = [value - mean_y for value in y]
    denominator = math.sqrt(sum(value * value for value in centered_x) * sum(value * value for value in centered_y))
    if denominator == 0.0:
        return None
    return sum(a * b for a, b in zip(centered_x, centered_y)) / denominator


def spearman(x: Sequence[float], y: Sequence[float]) -> float | None:
    return pearson(average_ranks(x), average_ranks(y))


def dominates(a: Sequence[float], b: Sequence[float]) -> bool:
    return all(left >= right for left, right in zip(a, b)) and any(
        left > right for left, right in zip(a, b)
    )


def pareto_fronts(vectors: Sequence[Sequence[float]]) -> list[int]:
    """Assign deterministic non-dominated fronts in O(n^2)."""
    n = len(vectors)
    domination_count = [0] * n
    dominated: list[list[int]] = [[] for _ in range(n)]
    for left in range(n):
        for right in range(left + 1, n):
            if dominates(vectors[left], vectors[right]):
                dominated[left].append(right)
                domination_count[right] += 1
            elif dominates(vectors[right], vectors[left]):
                dominated[right].append(left)
                domination_count[left] += 1

    assigned = [0] * n
    current = [index for index, count in enumerate(domination_count) if count == 0]
    front_number = 1
    while current:
        next_front: list[int] = []
        for index in current:
            assigned[index] = front_number
            for other in dominated[index]:
                domination_count[other] -= 1
                if domination_count[other] == 0:
                    next_front.append(other)
        current = next_front
        front_number += 1
    if any(front == 0 for front in assigned):
        raise RuntimeError("Pareto sorting did not assign every candidate")
    return assigned


def _index_rows(rows: Sequence[dict], label: str) -> dict[str, dict]:
    indexed: dict[str, dict] = {}
    for row in rows:
        candidate_id = row["candidate_id"]
        if candidate_id in indexed:
            raise ValueError(f"Duplicate candidate_id in {label}: {candidate_id}")
        indexed[candidate_id] = row
    return indexed


def _finite_probability(row: dict, column: str, label: str) -> float:
    try:
        value = float(row[column])
    except (TypeError, ValueError) as error:
        raise ValueError(f"Non-numeric {column} for {label}") from error
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"Invalid {column} for {label}: {value}")
    return value


def transform(
    layer1_rows: Sequence[dict],
    layer2_rows: Sequence[dict],
    layer3_rows: Sequence[dict],
    *,
    enforce_fixture: bool = False,
) -> dict[str, list[dict]]:
    require_columns(layer1_rows, IDENTITY_COLUMNS + ["layer1_accumulation_percentile"], "Layer 1")
    require_columns(layer2_rows, IDENTITY_COLUMNS + ["layer2_reference_score"], "Layer 2")
    require_columns(layer3_rows, IDENTITY_COLUMNS + ["layer3_reference_score"], "Layer 3")

    indexed1 = _index_rows(layer1_rows, "Layer 1")
    indexed2 = _index_rows(layer2_rows, "Layer 2")
    indexed3 = _index_rows(layer3_rows, "Layer 3")
    id_set = set(indexed1)
    if set(indexed2) != id_set or set(indexed3) != id_set:
        raise ValueError("Stage 09 candidate ID sets are not identical")

    ordered_ids = [row["candidate_id"] for row in layer1_rows]
    joined: list[dict] = []
    for candidate_id in ordered_ids:
        rows = (indexed1[candidate_id], indexed2[candidate_id], indexed3[candidate_id])
        for column in IDENTITY_COLUMNS:
            if len({row[column] for row in rows}) != 1:
                raise ValueError(f"Stage 09 identity mismatch for {candidate_id}: {column}")
        l1 = _finite_probability(rows[0], "layer1_accumulation_percentile", candidate_id)
        l2_score = _finite_probability(rows[1], "layer2_reference_score", candidate_id)
        l3_score = _finite_probability(rows[2], "layer3_reference_score", candidate_id)
        joined.append(
            {
                **{column: rows[0][column] for column in IDENTITY_COLUMNS},
                "layer1_accumulation_percentile": l1,
                "layer2_reference_score": l2_score,
                "layer3_reference_score": l3_score,
            }
        )

    counts = Counter(int(row["candidate_length_nt"]) for row in joined)
    if enforce_fixture and (len(joined) != 1375 or counts != {23: 688, 24: 687}):
        raise ValueError(f"Current Stage 10A accounting fixture failed: total={len(joined)}, lengths={dict(counts)}")

    groups: dict[tuple[str, int], list[int]] = defaultdict(list)
    for index, row in enumerate(joined):
        groups[(row["target_id"], int(row["candidate_length_nt"]))].append(index)

    correlations: list[dict] = []
    pareto_summary: list[dict] = []
    for (target_id, length), indexes in sorted(groups.items()):
        layer1 = [joined[index]["layer1_accumulation_percentile"] for index in indexes]
        layer2 = favourable_percentiles([joined[index]["layer2_reference_score"] for index in indexes])
        layer3 = favourable_percentiles([joined[index]["layer3_reference_score"] for index in indexes])
        equal_scores = [(a + b + c) / 3.0 for a, b, c in zip(layer1, layer2, layer3)]
        equal_ranks = average_ranks(equal_scores, descending=True)
        equal_percentiles = favourable_percentiles(equal_scores)
        fronts = pareto_fronts(list(zip(layer1, layer2, layer3)))

        for local_index, global_index in enumerate(indexes):
            row = joined[global_index]
            row.update(
                {
                    "stage10_layer1_percentile": layer1[local_index],
                    "stage10_layer2_percentile": layer2[local_index],
                    "stage10_layer3_percentile": layer3[local_index],
                    "stage10_equal_layer_score": equal_scores[local_index],
                    "stage10_equal_layer_rank": equal_ranks[local_index],
                    "stage10_equal_layer_percentile": equal_percentiles[local_index],
                    "stage10_pareto_front": fronts[local_index],
                    "stage10_minimum_layer_score": min(layer1[local_index], layer2[local_index], layer3[local_index]),
                }
            )

        for comparison, left_name, right_name, left, right in [
            ("L1_vs_L2", "stage10_layer1_percentile", "stage10_layer2_percentile", layer1, layer2),
            ("L1_vs_L3", "stage10_layer1_percentile", "stage10_layer3_percentile", layer1, layer3),
            ("L2_vs_L3", "stage10_layer2_percentile", "stage10_layer3_percentile", layer2, layer3),
        ]:
            rho = spearman(left, right)
            correlations.append(
                {
                    "target_id": target_id,
                    "candidate_length_nt": length,
                    "comparison": comparison,
                    "metric_x": left_name,
                    "metric_y": right_name,
                    "n_candidates": len(indexes),
                    "spearman_rho": "NA" if rho is None else rho,
                }
            )

        front_counts = Counter(fronts)
        for front in sorted(front_counts):
            pareto_summary.append(
                {
                    "target_id": target_id,
                    "candidate_length_nt": length,
                    "stage10_pareto_front": front,
                    "n_candidates": front_counts[front],
                    "fraction_candidates": front_counts[front] / len(indexes),
                }
            )

    return {"candidates": joined, "correlations": correlations, "pareto_summary": pareto_summary}


def qc_rows(result: dict[str, list[dict]], source_ids: set[str]) -> list[dict]:
    candidates = result["candidates"]
    counts = Counter(int(row["candidate_length_nt"]) for row in candidates)
    numeric = [
        "layer1_accumulation_percentile",
        "layer2_reference_score",
        "layer3_reference_score",
        "stage10_layer1_percentile",
        "stage10_layer2_percentile",
        "stage10_layer3_percentile",
        "stage10_equal_layer_score",
        "stage10_equal_layer_rank",
        "stage10_equal_layer_percentile",
        "stage10_pareto_front",
        "stage10_minimum_layer_score",
    ]
    finite = all(math.isfinite(float(row[column])) for row in candidates for column in numeric)
    arithmetic = max(
        abs(
            row["stage10_equal_layer_score"]
            - (row["stage10_layer1_percentile"] + row["stage10_layer2_percentile"] + row["stage10_layer3_percentile"]) / 3.0
        )
        for row in candidates
    )
    minimum = max(
        abs(
            row["stage10_minimum_layer_score"]
            - min(row["stage10_layer1_percentile"], row["stage10_layer2_percentile"], row["stage10_layer3_percentile"])
        )
        for row in candidates
    )
    columns = set(candidates[0])
    checks = [
        ("candidate_total", len(candidates) == 1375, len(candidates), "current Vd-CHIBIN fixture"),
        ("candidate_count_23nt", counts.get(23) == 688, counts.get(23, 0), "current Vd-CHIBIN fixture"),
        ("candidate_count_24nt", counts.get(24) == 687, counts.get(24, 0), "current Vd-CHIBIN fixture"),
        ("candidate_ids_preserved", {row["candidate_id"] for row in candidates} == source_ids, len(source_ids), "exact one-to-one join"),
        ("required_numeric_finite", finite, int(finite), "no NA/Inf"),
        ("layer1_exact_copy", all(row["stage10_layer1_percentile"] == row["layer1_accumulation_percentile"] for row in candidates), 0, "no reranking"),
        ("equal_layer_arithmetic", arithmetic == 0.0, arithmetic, "weights 1/3,1/3,1/3"),
        ("minimum_layer_arithmetic", minimum == 0.0, minimum, "minimum of L1/L2/L3"),
        ("pareto_fronts_exhaustive", all(int(row["stage10_pareto_front"]) >= 1 for row in candidates), len(candidates), "one positive front per candidate"),
        ("no_stage08_raw_metric", not (columns & FORBIDDEN_STAGE08_RAW_COLUMNS), len(columns & FORBIDDEN_STAGE08_RAW_COLUMNS), "Stage 09 layers only"),
        ("no_filter_or_gate", len(candidates) == len(source_ids), len(candidates), "all candidates retained"),
        ("no_region_or_construct_logic", not ({"region_score", "selected_region", "construct_score"} & columns), 0, "Stage 11 excluded"),
        ("no_efficacy_probability", "efficacy_probability" not in columns, 0, "not an efficacy model"),
    ]
    return [
        {"check": name, "status": "PASS" if passed else "FAIL", "value": value, "detail": detail}
        for name, passed, value, detail in checks
    ]


def run(layer1_path: Path, layer2_path: Path, layer3_path: Path, output_root: Path) -> None:
    layer1 = read_tsv(layer1_path)
    layer2 = read_tsv(layer2_path)
    layer3 = read_tsv(layer3_path)
    result = transform(layer1, layer2, layer3, enforce_fixture=True)
    qc = qc_rows(result, {row["candidate_id"] for row in layer1})
    failures = [row for row in qc if row["status"] == "FAIL"]
    if failures:
        raise ValueError("Stage 10A QC failed: " + "; ".join(row["check"] for row in failures))

    write_tsv(output_root / "candidate_stage10a.tsv", result["candidates"], CANDIDATE_COLUMNS)
    write_tsv(
        output_root / "stage10a_layer_correlations.tsv",
        result["correlations"],
        ["target_id", "candidate_length_nt", "comparison", "metric_x", "metric_y", "n_candidates", "spearman_rho"],
    )
    write_tsv(
        output_root / "stage10a_pareto_summary.tsv",
        result["pareto_summary"],
        ["target_id", "candidate_length_nt", "stage10_pareto_front", "n_candidates", "fraction_candidates"],
    )
    parameters = [
        {"parameter": "stage", "value": "10A", "detail": "individual-window integration only"},
        {"parameter": "ranking_stratum", "value": "target_id x candidate_length_nt", "detail": "lengths never pooled"},
        {"parameter": "percentile_formula", "value": "(average_ascending_rank - 0.5) / n", "detail": "higher is favourable"},
        {"parameter": "layer1_weight", "value": "0.3333333333333333", "detail": "neutral equal-layer reference"},
        {"parameter": "layer2_weight", "value": "0.3333333333333333", "detail": "neutral equal-layer reference"},
        {"parameter": "layer3_weight", "value": "0.3333333333333333", "detail": "neutral equal-layer reference"},
        {"parameter": "pareto_inputs", "value": "L1,L2,L3", "detail": "separate diagnostic; full precision"},
        {"parameter": "layer1_input", "value": str(layer1_path), "detail": "read-only Stage 09A output"},
        {"parameter": "layer2_input", "value": str(layer2_path), "detail": "read-only Stage 09B output"},
        {"parameter": "layer3_input", "value": str(layer3_path), "detail": "read-only Stage 09C output"},
    ]
    write_tsv(output_root / "stage10_parameters.tsv", parameters, ["parameter", "value", "detail"])
    write_tsv(output_root / "stage10_qc.tsv", qc, ["check", "status", "value", "detail"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layer1", type=Path, required=True)
    parser.add_argument("--layer2", type=Path, required=True)
    parser.add_argument("--layer3", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args.layer1, args.layer2, args.layer3, args.output_root)


if __name__ == "__main__":
    main()
