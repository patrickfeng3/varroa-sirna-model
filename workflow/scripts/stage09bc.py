#!/usr/bin/env python3
"""Canonical Stage 09B/09C deterministic transformations of Stage 08."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


IDENTIFIER_COLUMNS = (
    "target_id", "transcript_id", "display_name", "organism", "candidate_id",
    "candidate_length_nt", "start_1based", "end_1based", "target_sequence_dna",
    "target_sequence_rna", "antisense_guide_sequence_rna", "annotation_status",
    "start_region", "end_region", "overlap_regions", "crosses_annotation_boundary",
)
LAYER2_RAW = (
    "guide_5p_terminal_dg_4bp", "passenger_5p_terminal_dg_4bp", "asymmetry_ddg_4bp",
    "guide_5p_terminal_dg_5bp", "passenger_5p_terminal_dg_5bp", "asymmetry_ddg_5bp",
    "guide_self_fold_mfe_kcal_mol", "guide_self_fold_structure",
)
LAYER3_RAW = (
    "target_whole_p_unpaired", "target_whole_p_unpaired_w100_l80",
    "target_whole_p_unpaired_w200_l150", "target_seed_g2_8_p_unpaired",
    "target_seed_g2_8_p_unpaired_w100_l80", "target_seed_g2_8_p_unpaired_w200_l150",
)
EXPECTED_LENGTH_COUNTS = {23: 688, 24: 687}
FORBIDDEN_OUTPUT_COLUMNS = {
    "layer1_accumulation_linear_predictor", "layer1_accumulation_percentile",
    "overall_score", "overall_rank", "stage09_score", "stage09_rank",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


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
        raise ValueError("cannot percentile-normalize an empty group")
    return [(rank - 0.5) / len(values) for rank in average_ranks(values)]


def spearman_rho(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    x, y = average_ranks(left), average_ranks(right)
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    denominator = math.sqrt(
        sum((a - mean_x) ** 2 for a in x) * sum((b - mean_y) ** 2 for b in y)
    )
    return numerator / denominator if denominator else None


def _finite_float(row: Mapping[str, str], column: str) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid Stage 08 value for {column}: {row.get(column)!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite Stage 08 value for {column}")
    return value


def transform(rows: Sequence[Mapping[str, str]]) -> dict[str, list[dict[str, object]]]:
    if not rows:
        raise ValueError("Stage 08 candidate table is empty")
    required = set(IDENTIFIER_COLUMNS + LAYER2_RAW + LAYER3_RAW)
    if not required.issubset(rows[0]):
        raise ValueError(f"Stage 08 candidate schema lacks: {sorted(required.difference(rows[0]))}")
    if set(rows[0]).intersection(FORBIDDEN_OUTPUT_COLUMNS):
        raise ValueError("Stage 09A/overall fields are forbidden Stage 09B/09C inputs")
    seen: set[str] = set()
    counts: dict[int, int] = defaultdict(int)
    groups: dict[tuple[str, int], list[int]] = defaultdict(list)
    layer2 = []
    layer3 = []
    for index, source in enumerate(rows):
        candidate_id = source["candidate_id"]
        if candidate_id in seen:
            raise ValueError(f"duplicate Stage 08 candidate_id: {candidate_id}")
        seen.add(candidate_id)
        length = int(source["candidate_length_nt"])
        if length not in EXPECTED_LENGTH_COUNTS:
            raise ValueError(f"unsupported canonical Stage 09 candidate length: {length}")
        counts[length] += 1
        groups[(source["target_id"], length)].append(index)
        layer2.append({column: source[column] for column in IDENTIFIER_COLUMNS + LAYER2_RAW})
        layer3.append({column: source[column] for column in IDENTIFIER_COLUMNS + LAYER3_RAW})
    if len(rows) != 1375 or dict(counts) != EXPECTED_LENGTH_COUNTS:
        raise ValueError(f"Stage 09B/09C row-accounting failure: total={len(rows)}, lengths={dict(counts)}")

    layer2_metrics = {
        "asymmetry_4bp_percentile": "asymmetry_ddg_4bp",
        "asymmetry_5bp_percentile": "asymmetry_ddg_5bp",
        "guide_self_fold_percentile": "guide_self_fold_mfe_kcal_mol",
    }
    layer3_metrics = {
        "whole_site_accessibility_percentile": "target_whole_p_unpaired",
        "seed_accessibility_percentile": "target_seed_g2_8_p_unpaired",
    }
    for indexes in groups.values():
        for output_column, source_column in layer2_metrics.items():
            percentiles = favourable_percentiles([_finite_float(rows[index], source_column) for index in indexes])
            for index, percentile in zip(indexes, percentiles):
                layer2[index][output_column] = percentile
        for output_column, source_column in layer3_metrics.items():
            percentiles = favourable_percentiles([_finite_float(rows[index], source_column) for index in indexes])
            for index, percentile in zip(indexes, percentiles):
                layer3[index][output_column] = percentile

    for row in layer2:
        row["layer2_asymmetry_percentile"] = row["asymmetry_4bp_percentile"]
        row["layer2_self_fold_percentile"] = row["guide_self_fold_percentile"]
        row["asymmetry_4bp_5bp_percentile_difference"] = (
            row["asymmetry_4bp_percentile"] - row["asymmetry_5bp_percentile"]
        )
        row["layer2_component_difference"] = (
            row["layer2_asymmetry_percentile"] - row["layer2_self_fold_percentile"]
        )
        row["layer2_reference_score"] = 0.5 * (
            row["layer2_asymmetry_percentile"] + row["layer2_self_fold_percentile"]
        )
    for row in layer3:
        row["layer3_whole_accessibility_percentile"] = row["whole_site_accessibility_percentile"]
        row["layer3_seed_accessibility_percentile"] = row["seed_accessibility_percentile"]
        row["layer3_reference_score"] = 0.5 * (
            row["layer3_whole_accessibility_percentile"] + row["layer3_seed_accessibility_percentile"]
        )

    layer2_correlations = []
    layer3_correlations = []
    for (target_id, length), indexes in sorted(groups.items()):
        for left, right, name in (
            ("asymmetry_ddg_4bp", "guide_self_fold_mfe_kcal_mol", "asymmetry_4bp_vs_self_fold"),
            ("asymmetry_ddg_4bp", "asymmetry_ddg_5bp", "asymmetry_4bp_vs_5bp"),
        ):
            layer2_correlations.append({
                "target_id": target_id, "candidate_length_nt": length, "comparison": name,
                "metric_x": left, "metric_y": right, "n_candidates": len(indexes),
                "spearman_rho": spearman_rho(
                    [_finite_float(rows[index], left) for index in indexes],
                    [_finite_float(rows[index], right) for index in indexes],
                ),
            })
        comparisons = (
            ("target_whole_p_unpaired", "target_seed_g2_8_p_unpaired", "canonical_whole_vs_seed"),
            ("target_whole_p_unpaired", "target_whole_p_unpaired_w100_l80", "whole_canonical_vs_w100_l80"),
            ("target_whole_p_unpaired", "target_whole_p_unpaired_w200_l150", "whole_canonical_vs_w200_l150"),
            ("target_seed_g2_8_p_unpaired", "target_seed_g2_8_p_unpaired_w100_l80", "seed_canonical_vs_w100_l80"),
            ("target_seed_g2_8_p_unpaired", "target_seed_g2_8_p_unpaired_w200_l150", "seed_canonical_vs_w200_l150"),
        )
        for left, right, name in comparisons:
            layer3_correlations.append({
                "target_id": target_id, "candidate_length_nt": length, "comparison": name,
                "metric_x": left, "metric_y": right, "n_candidates": len(indexes),
                "spearman_rho": spearman_rho(
                    [_finite_float(rows[index], left) for index in indexes],
                    [_finite_float(rows[index], right) for index in indexes],
                ),
            })

    layer2_sensitivity = {23: [], 24: []}
    layer3_sensitivity = {23: [], 24: []}
    weights = [index / 10.0 for index in range(11)]
    for row2, row3 in zip(layer2, layer3):
        length = int(row2["candidate_length_nt"])
        for alpha in weights:
            layer2_sensitivity[length].append({
                "target_id": row2["target_id"], "candidate_id": row2["candidate_id"],
                "candidate_length_nt": length, "alpha": alpha,
                "layer2_alpha_score": alpha * row2["layer2_asymmetry_percentile"]
                + (1 - alpha) * row2["layer2_self_fold_percentile"],
                "status": "sensitivity_only_not_tuned",
            })
        for gamma in weights:
            layer3_sensitivity[length].append({
                "target_id": row3["target_id"], "candidate_id": row3["candidate_id"],
                "candidate_length_nt": length, "gamma": gamma,
                "layer3_gamma_score": gamma * row3["layer3_whole_accessibility_percentile"]
                + (1 - gamma) * row3["layer3_seed_accessibility_percentile"],
                "status": "sensitivity_only_not_tuned",
            })
    return {
        "layer2": layer2, "layer3": layer3,
        "layer2_correlations": layer2_correlations, "layer3_correlations": layer3_correlations,
        "layer2_sensitivity_23": layer2_sensitivity[23], "layer2_sensitivity_24": layer2_sensitivity[24],
        "layer3_sensitivity_23": layer3_sensitivity[23], "layer3_sensitivity_24": layer3_sensitivity[24],
    }


def run(stage08_candidates: Path, output_root: Path) -> None:
    result = transform(read_tsv(stage08_candidates))
    layer2_root = output_root / "09B_layer2_guide_competence"
    layer3_root = output_root / "09C_layer3_target_engagement"
    outputs = (
        (layer2_root / "candidate_layer2.tsv", result["layer2"]),
        (layer2_root / "layer2_weight_sensitivity_23nt.tsv", result["layer2_sensitivity_23"]),
        (layer2_root / "layer2_weight_sensitivity_24nt.tsv", result["layer2_sensitivity_24"]),
        (layer2_root / "layer2_correlations.tsv", result["layer2_correlations"]),
        (layer3_root / "candidate_layer3.tsv", result["layer3"]),
        (layer3_root / "layer3_weight_sensitivity_23nt.tsv", result["layer3_sensitivity_23"]),
        (layer3_root / "layer3_weight_sensitivity_24nt.tsv", result["layer3_sensitivity_24"]),
        (layer3_root / "layer3_correlations.tsv", result["layer3_correlations"]),
    )
    for path, rows in outputs:
        if not rows:
            raise ValueError(f"refusing to write empty Stage 09B/09C output: {path}")
        write_tsv(path, rows, list(rows[0]))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage08-candidates", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.stage08_candidates.is_file():
        raise FileNotFoundError(f"missing validated Stage 08 candidate table: {args.stage08_candidates}")
    run(args.stage08_candidates, args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
