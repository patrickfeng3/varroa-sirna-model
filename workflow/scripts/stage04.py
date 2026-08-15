#!/usr/bin/env python3
"""Canonical Stage 04 duplex-geometry aggregation and conditioned sequence analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import random
import statistics
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


STRANDS = ("sense", "antisense")
WEIGHTING_MODES = ("abundance", "unique_sequence")
POSITIONS = ("5p1", "5p2", "3p2", "3p1")
NUCLEOTIDES = ("A", "C", "G", "T")
JOINT_METRICS = (
    "varroa_2nt_joint_duplex_fraction",
    "varroa_2nt_reference_fraction_all",
    "varroa_2nt_reference_fraction_recovered",
    "varroa_2nt_reference_fraction_abundance_all",
    "varroa_2nt_reference_fraction_abundance_recovered",
)
RECOVERY_METRICS = (
    "passenger_recovery_fraction_unique",
    "passenger_recovery_fraction_abundance",
)
SEQUENCE_METRICS = (
    "joint_observed_fraction",
    "recovered_observed_fraction",
    "E_joint_absolute",
    "E_recovered_absolute",
    "E_all",
    "joint_vs_all_log2_contrast",
    "joint_vs_recovered_log2_contrast",
)


class Stage04Error(RuntimeError):
    """Structured Stage 04 input or consistency failure."""


@dataclass(frozen=True)
class Stage04Config:
    bootstrap_replicates: int
    random_seed: int
    ci_method: str
    ci_level: float
    frequency_sum_tolerance: float
    joint_5p_distance: int
    joint_3p_distance: int


def load_config(path: Path) -> Stage04Config:
    data = json.loads(path.read_text())["stage04"]
    config = Stage04Config(
        bootstrap_replicates=int(data["bootstrap_replicates"]),
        random_seed=int(data["random_seed"]),
        ci_method=str(data["ci_method"]),
        ci_level=float(data["ci_level"]),
        frequency_sum_tolerance=float(data["frequency_sum_tolerance"]),
        joint_5p_distance=int(data["joint_5p_distance"]),
        joint_3p_distance=int(data["joint_3p_distance"]),
    )
    if config.bootstrap_replicates <= 0:
        raise ValueError("Stage 04 bootstrap replicates must be positive")
    if config.ci_method != "percentile" or not 0 < config.ci_level < 1:
        raise ValueError("Stage 04 requires a percentile CI with 0 < level < 1")
    if (config.joint_5p_distance, config.joint_3p_distance) != (2, -2):
        raise ValueError("Stage 04 joint geometry must remain (+2,-2)")
    return config


def finite_value(value: object) -> float | None:
    if value is None or str(value) in {"", "NA", "NaN", "nan"}:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def median_or_none(values: Iterable[object]) -> float | None:
    usable = [x for value in values if (x := finite_value(value)) is not None]
    return statistics.median(usable) if usable else None


def safe_fraction(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def safe_log2_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or numerator <= 0 or denominator <= 0:
        return None
    return math.log2(numerator / denominator)


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def clustered_bootstrap(
    sample_values: dict[str, float],
    replicates: int,
    seed: int,
    level: float,
) -> tuple[float | None, float | None, int]:
    samples = sorted(sample_values)
    if not samples:
        return None, None, 0
    rng = random.Random(seed)
    estimates = [
        statistics.median(sample_values[rng.choice(samples)] for _ in samples)
        for _ in range(replicates)
    ]
    alpha = (1 - level) / 2
    return percentile(estimates, alpha), percentile(estimates, 1 - alpha), len(estimates)


def terminal_bases(sequence: str) -> dict[str, str]:
    if len(sequence) < 2:
        raise Stage04Error("focal sequence shorter than two nucleotides")
    return {"5p1": sequence[0], "5p2": sequence[1], "3p2": sequence[-2], "3p1": sequence[-1]}


def tied_ranks(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda x: x[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        rank = (index + 1 + end) / 2
        for original, _ in ordered[index:end]:
            ranks[original] = rank
        index = end
    return ranks


def spearman_rho(x: list[float], y: list[float]) -> float | None:
    if len(x) != len(y) or len(x) < 2:
        return None
    rx, ry = tied_ranks(x), tied_ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    numerator = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denominator = math.sqrt(
        sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)
    )
    return numerator / denominator if denominator else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def write_table(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def aggregate_metric_rows(
    rows: list[dict[str, object]],
    group_fields: tuple[str, ...],
    value_field: str,
    config: Stage04Config,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row[field] for field in group_fields)].append(row)
    sample_rows: list[dict[str, object]] = []
    across_rows: list[dict[str, object]] = []
    for key in sorted(grouped):
        group = grouped[key]
        by_sample: dict[str, list[float]] = defaultdict(list)
        undefined = 0
        for row in group:
            value = finite_value(row.get(value_field))
            if value is None:
                undefined += 1
            else:
                by_sample[str(row["sample"])].append(value)
        sample_values = {
            sample: statistics.median(values) for sample, values in by_sample.items() if values
        }
        common = dict(zip(group_fields, key))
        for sample in sorted(sample_values):
            sample_rows.append({
                **common,
                "sample": sample,
                "sample_median": sample_values[sample],
                "n_sample_virus_units": len(by_sample[sample]),
            })
        ci_low, ci_high, valid = clustered_bootstrap(
            sample_values, config.bootstrap_replicates, config.random_seed, config.ci_level
        )
        pair_values = [
            value for row in group if (value := finite_value(row.get(value_field))) is not None
        ]
        across_rows.append({
            **common,
            "sample_balanced_median": median_or_none(sample_values.values()),
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_samples": len(sample_values),
            "n_sample_virus_units": len(pair_values),
            "n_undefined_pair_values": undefined,
            "pair_balanced_median": median_or_none(pair_values),
            "bootstrap_replicates_requested": config.bootstrap_replicates,
            "bootstrap_replicates_valid": valid,
            "bootstrap_seed": config.random_seed,
            "ci_method": config.ci_method,
            "ci_level": config.ci_level,
        })
    return sample_rows, across_rows


def aggregate_full_spectrum(
    spectrum_rows: list[dict[str, str]], config: Stage04Config
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    normalized = []
    for row in spectrum_rows:
        for view, log_field, z_field in (
            ("duplex", "official_duplex_log_ratio", "official_duplex_wald_z"),
            (
                "unique_reference",
                "official_unique_reference_log_ratio",
                "official_unique_reference_wald_z",
            ),
        ):
            normalized.append({
                "sample": row["sample"],
                "analysis_unit": row["analysis_unit"],
                "biological_virus": row["biological_virus"],
                "focal_length": int(row["focal_length"]),
                "focal_strand": row["focal_strand"],
                "end": row["end"],
                "signed_distance": int(row["signed_distance"]),
                "official_view": view,
                "steprna_log_ratio": finite_value(row[log_field]),
                "steprna_wald_z": finite_value(row[z_field]),
            })
    groups = (
        "focal_length", "focal_strand", "end", "signed_distance", "official_view"
    )
    sample_rows, across = aggregate_metric_rows(normalized, groups, "steprna_log_ratio", config)
    z_sample, z_across = aggregate_metric_rows(normalized, groups, "steprna_wald_z", config)
    z_sample_lookup = {
        tuple(row[x] for x in (*groups, "sample")): row["sample_median"] for row in z_sample
    }
    z_across_lookup = {
        tuple(row[x] for x in groups): row["sample_balanced_median"] for row in z_across
    }
    for row in sample_rows:
        row["sample_steprna_log_ratio_median"] = row.pop("sample_median")
        row["sample_steprna_wald_z_median_descriptive"] = z_sample_lookup.get(
            tuple(row[x] for x in (*groups, "sample"))
        )
    for row in across:
        row["sample_balanced_steprna_log_ratio"] = row.pop("sample_balanced_median")
        row["sample_balanced_steprna_wald_z_descriptive"] = z_across_lookup.get(
            tuple(row[x] for x in groups)
        )
    strongest: dict[tuple[object, ...], float] = {}
    for row in across:
        value = row["sample_balanced_steprna_log_ratio"]
        if value is None:
            continue
        key = tuple(row[x] for x in ("focal_length", "focal_strand", "end", "official_view"))
        strongest[key] = max(strongest.get(key, -math.inf), float(value))
    for row in across:
        key = tuple(row[x] for x in ("focal_length", "focal_strand", "end", "official_view"))
        row["strongest_distance_indicator"] = int(
            row["sample_balanced_steprna_log_ratio"] is not None
            and row["sample_balanced_steprna_log_ratio"] == strongest.get(key)
        )
    return normalized, sample_rows, across


def subset_terminal_counts(
    focal_rows: list[dict[str, object]],
    subset_ids: set[str],
    weighting_mode: str,
) -> tuple[dict[tuple[str, str], float], float]:
    counts: dict[tuple[str, str], float] = defaultdict(float)
    total = 0.0
    for row in focal_rows:
        if str(row["focal_id"]) not in subset_ids:
            continue
        weight = float(row["focal_abundance"]) if weighting_mode == "abundance" else 1.0
        total += weight
        for position, nucleotide in terminal_bases(str(row["sequence"])).items():
            counts[(position, nucleotide)] += weight
    return counts, total


def calculate_sequence_pair_rows(
    focal_manifest: list[dict[str, str]],
    joint_manifest: list[dict[str, str]],
    recovered_by_run: dict[str, set[str]],
    expected_rows: list[dict[str, str]],
    general_rows: list[dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    focals_by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    focal_by_id: dict[str, dict[str, object]] = {}
    for row in focal_manifest:
        parsed = {**row, "focal_length": int(row["focal_length"]), "focal_abundance": float(row["focal_abundance"])}
        focal_id = row["focal_id"]
        if focal_id in focal_by_id:
            raise Stage04Error(f"duplicate focal ID in manifest: {focal_id}")
        focal_by_id[focal_id] = parsed
        focals_by_run[row["run_id"]].append(parsed)
    joint_by_run: dict[str, set[str]] = defaultdict(set)
    abundance_mismatches = 0
    absent_joint = 0
    for row in joint_manifest:
        focal_id = row["focal_id"]
        if focal_id not in focal_by_id:
            absent_joint += 1
            continue
        if not math.isclose(float(row["focal_abundance"]), float(focal_by_id[focal_id]["focal_abundance"])):
            abundance_mismatches += 1
        joint_by_run[row["run_id"]].add(focal_id)
    if absent_joint or abundance_mismatches:
        raise Stage04Error(
            f"joint/focal manifest inconsistency: absent={absent_joint}; abundance={abundance_mismatches}"
        )
    unknown_recovered = 0
    for run_id, recovered_ids in recovered_by_run.items():
        focal_ids = {str(row["focal_id"]) for row in focals_by_run.get(run_id, [])}
        unknown_recovered += len(recovered_ids - focal_ids)
    if unknown_recovered:
        raise Stage04Error(
            f"passenger-recovered references absent from focal manifest: {unknown_recovered}"
        )
    subset_violations = sum(
        len(joint_by_run[run_id] - recovered_by_run.get(run_id, set()))
        for run_id in joint_by_run
    )
    if subset_violations:
        raise Stage04Error(
            f"joint-support references outside recovered subset: {subset_violations}"
        )

    expected = {
        (
            row["sample"], row["analysis_unit"], int(row["length"]), row["strand_scope"],
            row["weighting_mode"], row["terminal_position"], row["nucleotide"],
        ): finite_value(row["expected_fraction"])
        for row in expected_rows if row["strand_scope"] in STRANDS
    }
    general = {
        (
            row["sample"], row["analysis_unit"], int(row["length"]), row["strand_scope"],
            row["weighting_mode"], row["terminal_position"], row["nucleotide"],
        ): finite_value(row["enrichment_ratio"])
        for row in general_rows if row["strand_scope"] in STRANDS
    }
    pair_rows = []
    missing_expected = missing_general = 0
    empty_joint = empty_recovered = 0
    max_frequency_deviation = 0.0
    for run_id in sorted(focals_by_run):
        focal_rows = focals_by_run[run_id]
        first = focal_rows[0]
        recovered_ids = recovered_by_run.get(run_id, set())
        joint_ids = joint_by_run.get(run_id, set())
        empty_joint += int(not joint_ids)
        empty_recovered += int(not recovered_ids)
        for mode in WEIGHTING_MODES:
            joint_counts, joint_total = subset_terminal_counts(focal_rows, joint_ids, mode)
            recovered_counts, recovered_total = subset_terminal_counts(focal_rows, recovered_ids, mode)
            for counts, total in ((joint_counts, joint_total), (recovered_counts, recovered_total)):
                if total:
                    for position in POSITIONS:
                        deviation = abs(sum(counts[(position, nt)] / total for nt in NUCLEOTIDES) - 1)
                        max_frequency_deviation = max(max_frequency_deviation, deviation)
            for position in POSITIONS:
                for nucleotide in NUCLEOTIDES:
                    key = (
                        first["sample"], first["analysis_unit"], int(first["focal_length"]),
                        first["focal_strand"], mode, position, nucleotide,
                    )
                    expected_fraction = expected.get(key)
                    all_enrichment = general.get(key)
                    missing_expected += int(key not in expected)
                    missing_general += int(key not in general)
                    joint_fraction = safe_fraction(joint_counts[(position, nucleotide)], joint_total)
                    recovered_fraction = safe_fraction(
                        recovered_counts[(position, nucleotide)], recovered_total
                    )
                    e_joint = (
                        joint_fraction / expected_fraction
                        if joint_fraction is not None and expected_fraction is not None
                        and expected_fraction > 0 else None
                    )
                    e_recovered = (
                        recovered_fraction / expected_fraction
                        if recovered_fraction is not None and expected_fraction is not None
                        and expected_fraction > 0 else None
                    )
                    pair_rows.append({
                        "sample": first["sample"],
                        "analysis_unit": first["analysis_unit"],
                        "biological_virus": first["biological_virus"],
                        "focal_length": int(first["focal_length"]),
                        "focal_strand": first["focal_strand"],
                        "weighting_mode": mode,
                        "terminal_position": position,
                        "nucleotide": nucleotide,
                        "n_focal_references": len(focal_rows),
                        "n_recovered_focal_references": len(recovered_ids),
                        "n_joint_focal_references": len(joint_ids),
                        "joint_total_weight": joint_total,
                        "recovered_total_weight": recovered_total,
                        "joint_observed_fraction": joint_fraction,
                        "recovered_observed_fraction": recovered_fraction,
                        "stage02_expected_fraction": expected_fraction,
                        "E_joint_absolute": e_joint,
                        "E_recovered_absolute": e_recovered,
                        "E_all": all_enrichment,
                        "joint_vs_all_log2_contrast": safe_log2_ratio(e_joint, all_enrichment),
                        "joint_vs_recovered_log2_contrast": safe_log2_ratio(e_joint, e_recovered),
                        "run_id": run_id,
                    })
    return pair_rows, {
        "missing_expected": missing_expected,
        "missing_general": missing_general,
        "empty_joint": empty_joint,
        "empty_recovered": empty_recovered,
        "max_frequency_deviation": max_frequency_deviation,
        "subset_violations": subset_violations,
        "abundance_mismatches": abundance_mismatches,
        "absent_joint": absent_joint,
    }


def aggregate_sequence_rows(
    pair_rows: list[dict[str, object]], config: Stage04Config
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    normalized = []
    keys = (
        "focal_length", "focal_strand", "weighting_mode", "terminal_position", "nucleotide"
    )
    for row in pair_rows:
        for metric in SEQUENCE_METRICS:
            normalized.append({
                "sample": row["sample"],
                "analysis_unit": row["analysis_unit"],
                **{field: row[field] for field in keys},
                "metric": metric,
                "value": row[metric],
            })
    sample_rows, across = aggregate_metric_rows(normalized, (*keys, "metric"), "value", config)
    for row in sample_rows:
        row["sample_metric_median"] = row.pop("sample_median")
    for row in across:
        row["sample_balanced_median"] = row["sample_balanced_median"]
    contrasts = [
        row for row in across
        if row["metric"] in {
            "joint_vs_all_log2_contrast", "joint_vs_recovered_log2_contrast"
        }
    ]
    return sample_rows, across, contrasts


def paired_comparisons(
    recovery_rows: list[dict[str, str]],
    joint_rows: list[dict[str, str]],
    normalized_spectrum: list[dict[str, object]],
    config: Stage04Config,
) -> list[dict[str, object]]:
    metric_rows = []
    for row in recovery_rows:
        for metric in RECOVERY_METRICS:
            metric_rows.append({
                "sample": row["sample"], "analysis_unit": row["analysis_unit"],
                "focal_length": int(row["focal_length"]), "focal_strand": row["focal_strand"],
                "metric": metric, "value": finite_value(row[metric]),
            })
    for row in joint_rows:
        for metric in JOINT_METRICS:
            metric_rows.append({
                "sample": row["sample"], "analysis_unit": row["analysis_unit"],
                "focal_length": int(row["focal_length"]), "focal_strand": row["focal_strand"],
                "metric": metric, "value": finite_value(row[metric]),
            })
    for row in normalized_spectrum:
        distance = int(row["signed_distance"])
        end = str(row["end"])
        if (end, distance) not in {("5p", 0), ("5p", 2), ("3p", 0), ("3p", -2)}:
            continue
        metric_rows.append({
            "sample": row["sample"], "analysis_unit": row["analysis_unit"],
            "focal_length": int(row["focal_length"]), "focal_strand": row["focal_strand"],
            "metric": f"steprna_{row['official_view']}_{end}_{distance:+d}_log_ratio",
            "value": row["steprna_log_ratio"],
        })
    values = {
        (row["sample"], row["analysis_unit"], row["focal_strand"], row["metric"], row["focal_length"]): row["value"]
        for row in metric_rows
    }
    pair_deltas = []
    prefixes = sorted({key[:4] for key in values})
    for sample, unit, strand, metric in prefixes:
        v23 = finite_value(values.get((sample, unit, strand, metric, 23)))
        v24 = finite_value(values.get((sample, unit, strand, metric, 24)))
        pair_deltas.append({
            "sample": sample, "analysis_unit": unit, "focal_strand": strand,
            "metric": metric,
            "paired_delta_24_minus_23": v24 - v23 if v23 is not None and v24 is not None else None,
        })
    _, across = aggregate_metric_rows(
        pair_deltas, ("focal_strand", "metric"), "paired_delta_24_minus_23", config
    )
    for row in across:
        row["sample_balanced_paired_delta_24_minus_23"] = row.pop("sample_balanced_median")
    return across


def calculate_redundancy(
    sequence_across: list[dict[str, object]],
    stage02_across: list[dict[str, str]],
) -> list[dict[str, object]]:
    lookup = {
        (
            int(row["focal_length"]), row["focal_strand"], row["weighting_mode"],
            row["terminal_position"], row["nucleotide"], row["metric"],
        ): finite_value(row["sample_balanced_median"])
        for row in sequence_across
    }
    general = {
        (
            int(row["length"]), row["strand_scope"], row["weighting_mode"],
            row["terminal_position"], row["nucleotide"],
        ): finite_value(row["sample_balanced_median_enrichment_ratio"])
        for row in stage02_across if row["strand_scope"] in STRANDS
    }
    output = []
    for length in (23, 24):
        for strand in STRANDS:
            for mode in WEIGHTING_MODES:
                pairs = []
                for position in POSITIONS:
                    for nucleotide in NUCLEOTIDES:
                        a = lookup.get((length, strand, mode, position, nucleotide, "E_joint_absolute"))
                        b = general.get((length, strand, mode, position, nucleotide))
                        if a is not None and b is not None:
                            pairs.append((a, b))
                output.append({
                    "focal_length": length, "focal_strand": strand,
                    "weighting_comparison": mode,
                    "comparison": "rho_joint_vs_general",
                    "n_matched_features": len(pairs),
                    "spearman_rho": spearman_rho(
                        [x[0] for x in pairs], [x[1] for x in pairs]
                    ),
                })
            pairs = []
            for position in POSITIONS:
                for nucleotide in NUCLEOTIDES:
                    a = lookup.get((
                        length, strand, "abundance", position, nucleotide,
                        "joint_vs_recovered_log2_contrast",
                    ))
                    b = lookup.get((
                        length, strand, "unique_sequence", position, nucleotide,
                        "joint_vs_recovered_log2_contrast",
                    ))
                    if a is not None and b is not None:
                        pairs.append((a, b))
            output.append({
                "focal_length": length, "focal_strand": strand,
                "weighting_comparison": "abundance_vs_unique_sequence",
                "comparison": "rho_joint_contrast_abundance_vs_unique",
                "n_matched_features": len(pairs),
                "spearman_rho": spearman_rho(
                    [x[0] for x in pairs], [x[1] for x in pairs]
                ),
            })
    return output


FIELDS = {
    "full_sample": [
        "sample", "focal_length", "focal_strand", "end", "signed_distance",
        "official_view", "sample_steprna_log_ratio_median",
        "sample_steprna_wald_z_median_descriptive", "n_sample_virus_units",
    ],
    "full_across": [
        "focal_length", "focal_strand", "end", "signed_distance", "official_view",
        "sample_balanced_steprna_log_ratio", "ci_low", "ci_high",
        "sample_balanced_steprna_wald_z_descriptive", "strongest_distance_indicator",
        "n_samples", "n_sample_virus_units", "n_undefined_pair_values",
        "pair_balanced_median", "bootstrap_replicates_requested",
        "bootstrap_replicates_valid", "bootstrap_seed", "ci_method", "ci_level",
    ],
    "metric_across": [
        "focal_length", "focal_strand", "metric", "sample_balanced_median",
        "ci_low", "ci_high", "n_samples", "n_sample_virus_units",
        "n_undefined_pair_values", "pair_balanced_median",
        "bootstrap_replicates_requested", "bootstrap_replicates_valid",
        "bootstrap_seed", "ci_method", "ci_level",
    ],
    "joint_sample": [
        "sample", "focal_length", "focal_strand", "metric", "sample_median",
        "n_sample_virus_units",
    ],
    "paired": [
        "focal_strand", "metric", "sample_balanced_paired_delta_24_minus_23",
        "ci_low", "ci_high", "n_samples", "n_sample_virus_units",
        "n_undefined_pair_values", "pair_balanced_median",
        "bootstrap_replicates_requested", "bootstrap_replicates_valid",
        "bootstrap_seed", "ci_method", "ci_level",
    ],
    "sequence_pair": [
        "sample", "analysis_unit", "biological_virus", "focal_length",
        "focal_strand", "weighting_mode", "terminal_position", "nucleotide",
        "n_focal_references", "n_recovered_focal_references",
        "n_joint_focal_references", "joint_total_weight", "recovered_total_weight",
        "joint_observed_fraction", "recovered_observed_fraction",
        "stage02_expected_fraction", "E_joint_absolute", "E_recovered_absolute",
        "E_all", "joint_vs_all_log2_contrast",
        "joint_vs_recovered_log2_contrast", "run_id",
    ],
    "sequence_sample": [
        "sample", "focal_length", "focal_strand", "weighting_mode",
        "terminal_position", "nucleotide", "metric", "sample_metric_median",
        "n_sample_virus_units",
    ],
    "sequence_across": [
        "focal_length", "focal_strand", "weighting_mode", "terminal_position",
        "nucleotide", "metric", "sample_balanced_median", "ci_low", "ci_high",
        "n_samples", "n_sample_virus_units", "n_undefined_pair_values",
        "pair_balanced_median", "bootstrap_replicates_requested",
        "bootstrap_replicates_valid", "bootstrap_seed", "ci_method", "ci_level",
    ],
    "redundancy": [
        "focal_length", "focal_strand", "weighting_comparison", "comparison",
        "n_matched_features", "spearman_rho",
    ],
    "qc": ["metric", "status", "value", "details"],
}


def run_stage04(
    stage02_root: Path, stage03_root: Path, config_path: Path, output_root: Path
) -> tuple[float, bool]:
    started = time.monotonic()
    config = load_config(config_path)
    recovery = read_tsv(stage03_root / "parsed/passenger_recovery_by_pair.tsv")
    spectrum = read_tsv(stage03_root / "parsed/overhang_spectrum_by_pair.tsv")
    joint = read_tsv(stage03_root / "parsed/joint_geometry_by_pair.tsv")
    joint_refs = read_gzip_tsv(stage03_root / "parsed/joint_geometry_references.tsv.gz")
    focals = read_gzip_tsv(stage03_root / "inputs/focal_reference_manifest.tsv.gz")
    run_manifest = read_tsv(stage03_root / "provenance/run_manifest.tsv")
    expected = read_tsv(stage02_root / "background/terminal_expected_by_pair.tsv")
    general_pair = read_tsv(stage02_root / "enrichment/terminal_enrichment_by_pair.tsv")
    general_across = read_tsv(stage02_root / "enrichment/terminal_enrichment_across_dataset.tsv")

    recovered_by_run: dict[str, set[str]] = {}
    raw_missing = 0
    for run in run_manifest:
        run_id = run["run_id"]
        path = stage03_root / "raw" / run_id / f"{run_id}_passenger_number.csv"
        if not path.exists():
            raw_missing += 1
            recovered_by_run[run_id] = set()
            continue
        with path.open(newline="") as handle:
            recovered_by_run[run_id] = {
                row["siRNA_reference"] for row in csv.DictReader(handle)
                if int(row["number"]) > 0
            }
    if raw_missing:
        raise Stage04Error(f"missing canonical Stage 03 passenger-number files: {raw_missing}")

    normalized_spectrum, full_sample, full_across = aggregate_full_spectrum(spectrum, config)
    recovery_normalized = [
        {
            "sample": row["sample"], "analysis_unit": row["analysis_unit"],
            "focal_length": int(row["focal_length"]), "focal_strand": row["focal_strand"],
            "metric": metric, "value": finite_value(row[metric]),
        }
        for row in recovery for metric in RECOVERY_METRICS
    ]
    _, recovery_across = aggregate_metric_rows(
        recovery_normalized, ("focal_length", "focal_strand", "metric"), "value", config
    )
    joint_normalized = [
        {
            "sample": row["sample"], "analysis_unit": row["analysis_unit"],
            "focal_length": int(row["focal_length"]), "focal_strand": row["focal_strand"],
            "metric": metric, "value": finite_value(row[metric]),
        }
        for row in joint for metric in JOINT_METRICS
    ]
    joint_sample, joint_across = aggregate_metric_rows(
        joint_normalized, ("focal_length", "focal_strand", "metric"), "value", config
    )
    paired = paired_comparisons(recovery, joint, normalized_spectrum, config)
    sequence_pair, sequence_qc = calculate_sequence_pair_rows(
        focals, joint_refs, recovered_by_run, expected, general_pair
    )
    sequence_sample, sequence_across, contrasts = aggregate_sequence_rows(sequence_pair, config)
    redundancy = calculate_redundancy(sequence_across, general_across)

    write_table(output_root / "population/full_spectrum_by_sample.tsv", full_sample, FIELDS["full_sample"])
    write_table(output_root / "population/full_spectrum_across_dataset.tsv", full_across, FIELDS["full_across"])
    write_table(output_root / "population/passenger_recovery_across_dataset.tsv", recovery_across, FIELDS["metric_across"])
    write_table(output_root / "population/joint_geometry_by_sample.tsv", joint_sample, FIELDS["joint_sample"])
    write_table(output_root / "population/joint_geometry_across_dataset.tsv", joint_across, FIELDS["metric_across"])
    write_table(output_root / "comparisons/paired_23_vs_24.tsv", paired, FIELDS["paired"])
    write_table(output_root / "sequence_features/geometry_terminal_by_pair.tsv", sequence_pair, FIELDS["sequence_pair"])
    write_table(output_root / "sequence_features/geometry_terminal_by_sample.tsv", sequence_sample, FIELDS["sequence_sample"])
    write_table(output_root / "sequence_features/geometry_terminal_across_dataset.tsv", sequence_across, FIELDS["sequence_across"])
    write_table(output_root / "sequence_features/geometry_specific_contrasts.tsv", contrasts, FIELDS["sequence_across"])
    write_table(output_root / "sequence_features/redundancy.tsv", redundancy, FIELDS["redundancy"])

    qc: list[dict[str, object]] = []
    def q(metric: str, value: object, status: str = "PASS", details: str = "") -> None:
        qc.append({"metric": metric, "status": status, "value": value, "details": details})
    run_keys = [
        (row["sample"], row["analysis_unit"], int(row["focal_length"]), row["focal_strand"])
        for row in run_manifest
    ]
    q("stage03_runs_represented", len(run_manifest))
    q("samples_represented", len({x[0] for x in run_keys}))
    q("sample_virus_units_represented", len({x[:2] for x in run_keys}))
    classes = sorted({(x[2], x[3]) for x in run_keys})
    q("focal_classes_represented", len(classes), "PASS" if len(classes) == 4 else "FAIL", str(classes))
    duplicate_keys = len(run_keys) - len(set(run_keys))
    q("duplicate_fixed_run_keys", duplicate_keys, "FAIL" if duplicate_keys else "PASS")
    expected_keys = len({x[:2] for x in run_keys}) * 4
    missing_keys = expected_keys - len(set(run_keys))
    q("missing_fixed_run_keys", missing_keys, "FAIL" if missing_keys else "PASS")
    required_spectrum = {
        (row["sample"], row["analysis_unit"], int(row["focal_length"]), row["focal_strand"], row["end"], int(row["signed_distance"]))
        for row in spectrum
    }
    missing_required = 0
    for sample, unit in sorted({x[:2] for x in run_keys}):
        for length in (23, 24):
            for strand in STRANDS:
                for end, distance in (("5p", 0), ("5p", 2), ("3p", 0), ("3p", -2)):
                    missing_required += int((sample, unit, length, strand, end, distance) not in required_spectrum)
    q("missing_required_signed_distance_rows", missing_required, "WARN" if missing_required else "PASS")
    undefined_logs = sum(row["steprna_log_ratio"] is None for row in normalized_spectrum)
    q("missing_or_undefined_official_log_ratios", undefined_logs, "INFO")
    q("joint_support_references_absent_from_focal_manifest", sequence_qc["absent_joint"], "FAIL" if sequence_qc["absent_joint"] else "PASS")
    q("joint_support_references_not_recovered", sequence_qc["subset_violations"], "FAIL" if sequence_qc["subset_violations"] else "PASS")
    q("focal_abundance_mismatches", sequence_qc["abundance_mismatches"], "FAIL" if sequence_qc["abundance_mismatches"] else "PASS")
    q("missing_stage02_expected_matches", sequence_qc["missing_expected"], "FAIL" if sequence_qc["missing_expected"] else "PASS")
    q("missing_stage02_general_enrichment_matches", sequence_qc["missing_general"], "FAIL" if sequence_qc["missing_general"] else "PASS")
    deviation = float(sequence_qc["max_frequency_deviation"])
    q("maximum_terminal_frequency_sum_deviation", deviation, "PASS" if deviation <= config.frequency_sum_tolerance else "FAIL")
    q("empty_joint_support_subsets", sequence_qc["empty_joint"], "WARN" if sequence_qc["empty_joint"] else "PASS")
    q("empty_passenger_recovered_subsets", sequence_qc["empty_recovered"], "WARN" if sequence_qc["empty_recovered"] else "PASS")
    undefined_contrasts = sum(
        row[metric] is None
        for row in sequence_pair
        for metric in ("joint_vs_all_log2_contrast", "joint_vs_recovered_log2_contrast")
    )
    q("undefined_pair_level_log2_contrasts", undefined_contrasts, "INFO")
    all_across = [*full_across, *recovery_across, *joint_across, *sequence_across, *paired]
    sample_counts = [int(row["n_samples"]) for row in all_across]
    q("canonical_summary_contributing_samples_min", min(sample_counts, default=0))
    q("canonical_summary_contributing_samples_max", max(sample_counts, default=0))
    valid_bootstraps = [int(row["bootstrap_replicates_valid"]) for row in all_across]
    q("bootstrap_replicates_requested", config.bootstrap_replicates)
    q("bootstrap_replicates_valid_min", min(valid_bootstraps, default=0))
    q("bootstrap_replicates_valid_max", max(valid_bootstraps, default=0))
    q("historical_delta_dicer_status", "deferred_not_run", "INFO")
    write_table(output_root / "qc/stage04_accounting.tsv", qc, FIELDS["qc"])
    return time.monotonic() - started, any(row["status"] == "FAIL" for row in qc)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage02-root", required=True, type=Path)
    parser.add_argument("--stage03-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        elapsed, failed = run_stage04(
            args.stage02_root.resolve(), args.stage03_root.resolve(),
            args.config.resolve(), args.output_root.resolve(),
        )
        print(f"Stage 04 completed in {elapsed:.3f} seconds", file=sys.stderr)
        return 1 if failed else 0
    except Exception as exc:
        output = args.output_root.resolve() / "qc/stage04_accounting.tsv"
        write_table(
            output,
            [{"metric": "stage04_execution", "status": "FAIL", "value": 1, "details": str(exc)}],
            FIELDS["qc"],
        )
        print(f"Stage 04 failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
