#!/usr/bin/env python3
"""Aggregate the canonical Stage 03 same-duplex joint geometry spectrum."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

from stage04 import Stage04Config, aggregate_metric_rows, load_config, write_table


PAIR_KEYS = ("sample", "analysis_unit", "biological_virus", "focal_length", "focal_strand", "run_id")
GEOMETRY_KEYS = ("steprna_5p_distance", "steprna_3p_distance")
GROUP_KEYS = ("focal_length", "focal_strand", *GEOMETRY_KEYS)
INTERPRETATION = (
    "Marginal end-distance 0 is prominent/recurrent; fully blunt (0,0) duplexes "
    "are a minority; pre-specified (+2,-2) duplexes are a minority; the joint "
    "duplex landscape is heterogeneous."
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_pair_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    return [
        {
            **{key: row[key] for key in PAIR_KEYS},
            "focal_length": int(row["focal_length"]),
            "steprna_5p_distance": int(row["steprna_5p_distance"]),
            "steprna_3p_distance": int(row["steprna_3p_distance"]),
            "official_duplex_count": int(row["official_duplex_count"]),
            "total_recovered_duplexes": int(row["total_recovered_duplexes"]),
            "joint_duplex_fraction": float(row["joint_duplex_fraction"]),
        }
        for row in rows
    ]


def densify_joint_spectrum(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Represent absent geometries as zero within every non-empty biological run."""
    by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    geometries: dict[tuple[int, str], set[tuple[int, int]]] = defaultdict(set)
    for row in rows:
        by_run[str(row["run_id"])].append(row)
        class_key = (int(row["focal_length"]), str(row["focal_strand"]))
        geometries[class_key].add(
            (int(row["steprna_5p_distance"]), int(row["steprna_3p_distance"]))
        )
    for class_key in geometries:
        geometries[class_key].update({(0, 0), (2, -2)})

    dense: list[dict[str, object]] = []
    for run_id in sorted(by_run):
        run_rows = by_run[run_id]
        first = run_rows[0]
        class_key = (int(first["focal_length"]), str(first["focal_strand"]))
        counts = {
            (int(row["steprna_5p_distance"]), int(row["steprna_3p_distance"])):
            int(row["official_duplex_count"])
            for row in run_rows
        }
        total = int(first["total_recovered_duplexes"])
        for d5, d3 in sorted(geometries[class_key]):
            count = counts.get((d5, d3), 0)
            dense.append({
                **{key: first[key] for key in PAIR_KEYS},
                "steprna_5p_distance": d5,
                "steprna_3p_distance": d3,
                "official_duplex_count": count,
                "total_recovered_duplexes": total,
                "joint_duplex_fraction": count / total,
            })
    return dense


def summarize_modes(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_run[str(row["run_id"])].append(row)
    output = []
    for run_id in sorted(by_run):
        group = by_run[run_id]
        maximum = max(int(row["official_duplex_count"]) for row in group)
        tied = sorted(
            (
                int(row["steprna_5p_distance"]),
                int(row["steprna_3p_distance"]),
            )
            for row in group
            if int(row["official_duplex_count"]) == maximum
        )
        zero_row = next(
            (row for row in group if int(row["steprna_5p_distance"]) == 0 and int(row["steprna_3p_distance"]) == 0),
            None,
        )
        first = group[0]
        output.append({
            **{key: first[key] for key in PAIR_KEYS},
            "most_common_steprna_5p_distance": tied[0][0],
            "most_common_steprna_3p_distance": tied[0][1],
            "most_common_official_duplex_count": maximum,
            "n_tied_most_common_geometries": len(tied),
            "zero_zero_is_most_common": int((0, 0) in tied),
            "zero_zero_joint_duplex_fraction": (
                float(zero_row["joint_duplex_fraction"]) if zero_row else 0.0
            ),
            "zero_zero_fraction_gt_0_5": int(
                zero_row is not None and float(zero_row["joint_duplex_fraction"]) > 0.5
            ),
        })
    return output


def aggregate_joint_spectrum(
    sparse_rows: list[dict[str, object]], config: Stage04Config
) -> tuple[
    list[dict[str, object]], list[dict[str, object]], list[dict[str, object]],
    list[dict[str, object]], list[dict[str, object]], list[dict[str, object]],
]:
    dense = densify_joint_spectrum(sparse_rows)
    sample_rows, across_rows = aggregate_metric_rows(
        dense, GROUP_KEYS, "joint_duplex_fraction", config
    )
    for row in sample_rows:
        row["sample_joint_duplex_fraction_median"] = row.pop("sample_median")

    dense_groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in dense:
        dense_groups[tuple(row[key] for key in GROUP_KEYS)].append(row)
    for row in across_rows:
        row["sample_balanced_median_joint_duplex_fraction"] = row.pop("sample_balanced_median")
        group = dense_groups[tuple(row[key] for key in GROUP_KEYS)]
        pooled_count = sum(int(item["official_duplex_count"]) for item in group)
        pooled_total = sum(int(item["total_recovered_duplexes"]) for item in group)
        row["pooled_official_duplex_count"] = pooled_count
        row["pooled_total_recovered_duplexes"] = pooled_total
        row["pooled_joint_duplex_fraction"] = pooled_count / pooled_total

    modes = summarize_modes(sparse_rows)
    lookup = {
        (int(row["focal_length"]), str(row["focal_strand"]), int(row["steprna_5p_distance"]), int(row["steprna_3p_distance"])): row
        for row in across_rows
    }
    summary = []
    represented_classes = sorted(
        {(int(row["focal_length"]), str(row["focal_strand"])) for row in sparse_rows}
    )
    for focal_length, focal_strand in represented_classes:
            zero = lookup[(focal_length, focal_strand, 0, 0)]
            varroa = lookup[(focal_length, focal_strand, 2, -2)]
            class_modes = [
                row for row in modes
                if int(row["focal_length"]) == focal_length and row["focal_strand"] == focal_strand
            ]
            summary.append({
                "focal_length": focal_length,
                "focal_strand": focal_strand,
                "sample_balanced_median_zero_zero_fraction": zero["sample_balanced_median_joint_duplex_fraction"],
                "zero_zero_ci_low": zero["ci_low"],
                "zero_zero_ci_high": zero["ci_high"],
                "pair_balanced_median_zero_zero_fraction": zero["pair_balanced_median"],
                "pooled_zero_zero_fraction": zero["pooled_joint_duplex_fraction"],
                "sample_balanced_median_plus2_minus2_fraction": varroa["sample_balanced_median_joint_duplex_fraction"],
                "plus2_minus2_ci_low": varroa["ci_low"],
                "plus2_minus2_ci_high": varroa["ci_high"],
                "pair_balanced_median_plus2_minus2_fraction": varroa["pair_balanced_median"],
                "pooled_plus2_minus2_fraction": varroa["pooled_joint_duplex_fraction"],
                "n_runs": len(class_modes),
                "n_runs_zero_zero_most_common": sum(int(row["zero_zero_is_most_common"]) for row in class_modes),
                "fraction_runs_zero_zero_most_common": sum(int(row["zero_zero_is_most_common"]) for row in class_modes) / len(class_modes),
                "n_runs_zero_zero_fraction_gt_0_5": sum(int(row["zero_zero_fraction_gt_0_5"]) for row in class_modes),
                "interpretation": INTERPRETATION,
            })

    by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in sparse_rows:
        by_run[str(row["run_id"])].append(row)
    sum_bad = 0
    max_fraction_deviation = 0.0
    for group in by_run.values():
        total = int(group[0]["total_recovered_duplexes"])
        count_sum = sum(int(row["official_duplex_count"]) for row in group)
        sum_bad += int(count_sum != total)
        max_fraction_deviation = max(
            max_fraction_deviation,
            abs(count_sum / total - 1.0),
        )
    qc = [
        {"metric": "joint_geometry_runs", "status": "PASS", "value": len(by_run), "details": ""},
        {"metric": "joint_geometry_count_sum_inconsistencies", "status": "FAIL" if sum_bad else "PASS", "value": sum_bad, "details": "counts must sum exactly to total_recovered_duplexes per run"},
        {"metric": "maximum_joint_fraction_sum_deviation", "status": "FAIL" if max_fraction_deviation > config.frequency_sum_tolerance else "PASS", "value": max_fraction_deviation, "details": f"tolerance={config.frequency_sum_tolerance}"},
    ]
    return dense, sample_rows, across_rows, modes, summary, qc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage03-spectrum", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args(argv)
    started = time.monotonic()
    config = load_config(args.config.resolve())
    sparse = normalize_pair_rows(read_tsv(args.stage03_spectrum.resolve()))
    _, samples, across, modes, summary, qc = aggregate_joint_spectrum(sparse, config)
    root = args.output_root.resolve()
    write_table(root / "population/joint_geometry_spectrum_by_sample.tsv", samples, [*GROUP_KEYS, "sample", "sample_joint_duplex_fraction_median", "n_sample_virus_units"])
    write_table(root / "population/joint_geometry_spectrum_across_dataset.tsv", across, [*GROUP_KEYS, "sample_balanced_median_joint_duplex_fraction", "ci_low", "ci_high", "n_samples", "n_sample_virus_units", "n_undefined_pair_values", "pair_balanced_median", "pooled_official_duplex_count", "pooled_total_recovered_duplexes", "pooled_joint_duplex_fraction", "bootstrap_replicates_requested", "bootstrap_replicates_valid", "bootstrap_seed", "ci_method", "ci_level"])
    write_table(root / "population/joint_geometry_mode_by_pair.tsv", modes, [*PAIR_KEYS, "most_common_steprna_5p_distance", "most_common_steprna_3p_distance", "most_common_official_duplex_count", "n_tied_most_common_geometries", "zero_zero_is_most_common", "zero_zero_joint_duplex_fraction", "zero_zero_fraction_gt_0_5"])
    write_table(root / "population/joint_geometry_spectrum_summary.tsv", summary, ["focal_length", "focal_strand", "sample_balanced_median_zero_zero_fraction", "zero_zero_ci_low", "zero_zero_ci_high", "pair_balanced_median_zero_zero_fraction", "pooled_zero_zero_fraction", "sample_balanced_median_plus2_minus2_fraction", "plus2_minus2_ci_low", "plus2_minus2_ci_high", "pair_balanced_median_plus2_minus2_fraction", "pooled_plus2_minus2_fraction", "n_runs", "n_runs_zero_zero_most_common", "fraction_runs_zero_zero_most_common", "n_runs_zero_zero_fraction_gt_0_5", "interpretation"])
    write_table(root / "qc/stage04_joint_geometry_spectrum_accounting.tsv", qc, ["metric", "status", "value", "details"])
    print(f"Stage 04 joint-spectrum aggregation completed in {time.monotonic()-started:.3f} seconds", file=sys.stderr)
    return 1 if any(row["status"] == "FAIL" for row in qc) else 0


if __name__ == "__main__":
    sys.exit(main())
