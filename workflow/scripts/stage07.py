#!/usr/bin/env python3
"""Canonical Stage 07 empirical guide-sequence association landscape."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import random
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Sequence


NUCLEOTIDES = ("A", "C", "G", "T")
STRANDS = ("sense", "antisense")
TARGET_LENGTHS = (23, 24)
REPRESENTATION_MODES = ("unique_sequence", "abundance")
ALL_MODES = REPRESENTATION_MODES + ("accumulation",)
BOOTSTRAP_REPLICATES = 5000
BOOTSTRAP_SEED = 20260814
CI_LEVEL = 0.95
FREQUENCY_TOLERANCE = 1e-12
REGRESSION_TOLERANCE = 1e-12
EXPECTED_PRIMARY_SAMPLES = 20
EXPECTED_PRIMARY_PAIRS = 54
REGIONAL_GC_WIDTH = 6
COMPLEMENT = str.maketrans("ACGT", "TGCA")
EXPECTED_CATEGORIES = {
    "mapping_mode": {"exact", "1mm"},
    "virus_assignment": {"assigned", "ambiguous_multi_virus"},
    "strand": {"sense", "antisense", "ambiguous"},
}


PAIR_FIELDS = [
    "sample", "analysis_unit", "biological_virus", "polarity", "length", "strand",
    "weighting_mode", "position_5p", "position_from_3p", "nucleotide",
    "observed_nucleotide_weight", "observed_total_weight", "observed_fraction",
    "expected_nucleotide_weight", "expected_total_windows", "expected_fraction",
    "representation_enrichment", "representation_delta_fraction",
    "unique_fraction", "abundance_fraction",
    "accumulation_ratio", "accumulation_delta_fraction", "log2_accumulation_ratio",
    "valid_background_windows",
]

SAMPLE_FIELDS = [
    "sample", "length", "strand", "weighting_mode", "position_5p", "position_from_3p",
    "nucleotide", "n_virus_units", "sample_observed_fraction", "sample_expected_fraction",
    "sample_representation_enrichment", "sample_representation_delta_fraction",
    "sample_unique_fraction", "sample_abundance_fraction", "sample_accumulation_ratio",
    "sample_accumulation_delta_fraction", "sample_log2_accumulation_ratio",
]

SUMMARY_FIELDS = [
    "strand", "length", "position_5p", "position_from_3p", "nucleotide",
    "weighting_mode", "endpoint", "n_samples", "n_samples_total",
    "sample_balanced_observed_fraction", "sample_balanced_expected_fraction",
    "sample_balanced_representation_enrichment",
    "sample_balanced_representation_delta_fraction", "sample_balanced_unique_fraction",
    "sample_balanced_abundance_fraction", "sample_balanced_accumulation_ratio",
    "sample_balanced_accumulation_delta_fraction", "log2_accumulation_ratio",
    "bootstrap_ci_low", "bootstrap_ci_high", "bootstrap_replicates_requested",
    "bootstrap_replicates_valid", "bootstrap_seed", "ci_method", "ci_level",
    "sign_test_n_nonzero", "sign_test_n_positive", "sign_test_n_negative",
    "sign_test_estimability", "raw_p", "bh_p", "by_p",
]

GC_PAIR_FIELDS = [
    "sample", "analysis_unit", "biological_virus", "polarity", "length", "strand",
    "observed_GC9_14_mean_unique", "observed_GC9_14_mean_abundance",
    "expected_GC9_14_mean", "GC9_14_delta_unique_vs_expected",
    "GC9_14_delta_abundance_vs_expected", "GC9_14_accumulation_delta",
]

GC_SAMPLE_FIELDS = [
    "sample", "length", "strand", "n_virus_units",
    "observed_GC9_14_mean_unique", "observed_GC9_14_mean_abundance",
    "expected_GC9_14_mean", "GC9_14_delta_unique_vs_expected",
    "GC9_14_delta_abundance_vs_expected", "GC9_14_accumulation_delta",
]

GC_SUMMARY_FIELDS = [
    "strand", "length", "endpoint", "effect_metric", "n_samples", "effect_estimate",
    "bootstrap_ci_low", "bootstrap_ci_high", "bootstrap_replicates_requested",
    "bootstrap_replicates_valid", "bootstrap_seed", "ci_method", "ci_level",
    "sign_test_n_nonzero", "sign_test_n_positive", "sign_test_n_negative",
    "sign_test_estimability", "raw_p", "bh_p",
]

LITERATURE_FIELDS = [
    "feature", "strand", "length", "endpoint", "effect_metric", "effect_estimate",
    "paired_delta_estimate", "bootstrap_ci_low", "bootstrap_ci_high", "n_samples",
    "sign_test_n_nonzero", "sign_test_n_positive", "sign_test_n_negative",
    "sign_test_estimability", "raw_p", "bh_p",
]

REGIONAL_PAIR_FIELDS = [
    "sample", "analysis_unit", "biological_virus", "polarity", "strand", "length",
    "start_5p", "end_5p", "near_3p", "far_3p", "region_5p", "region_3p",
    "observed_gc6_mean_unique", "observed_gc6_mean_abundance", "expected_gc6_mean",
    "regional_gc6_delta_unique_vs_expected",
    "regional_gc6_delta_abundance_vs_expected", "regional_gc6_accumulation_delta",
]

REGIONAL_SAMPLE_FIELDS = [
    "sample", "strand", "length", "start_5p", "end_5p", "near_3p", "far_3p",
    "region_5p", "region_3p", "n_virus_units", "observed_gc6_mean_unique",
    "observed_gc6_mean_abundance", "expected_gc6_mean",
    "regional_gc6_delta_unique_vs_expected",
    "regional_gc6_delta_abundance_vs_expected", "regional_gc6_accumulation_delta",
]

REGIONAL_SUMMARY_FIELDS = [
    "strand", "length", "start_5p", "end_5p", "near_3p", "far_3p",
    "region_5p", "region_3p", "endpoint", "n_samples",
    "sample_balanced_observed_gc6", "sample_balanced_expected_gc6",
    "sample_balanced_regional_gc6_delta", "bootstrap_ci_low", "bootstrap_ci_high",
    "bootstrap_replicates_requested", "bootstrap_replicates_valid", "bootstrap_seed",
    "ci_method", "ci_level", "sign_test_n_nonzero", "sign_test_n_positive",
    "sign_test_n_negative", "sign_test_estimability", "raw_p", "evidence_class",
    "validation_bh_p", "regional_bh_p", "regional_by_p",
]

QC_FIELDS = ["metric", "status", "value", "details"]
REGRESSION_FIELDS = [
    "level", "sample", "analysis_unit", "length", "strand", "weighting_mode",
    "terminal_position", "position_5p", "nucleotide", "metric", "stage02_value",
    "stage07_value", "absolute_difference", "status",
]


class Stage07Error(RuntimeError):
    """A canonical Stage 07 validation failure."""


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def observed_physical_sequence(sequence: str, strand: str) -> str:
    """Observed reads already use physical sequenced 5'-to-3' orientation."""
    if strand not in STRANDS:
        raise ValueError(f"invalid observed strand: {strand}")
    return sequence


def expected_physical_sequence(reference_window: str, strand: str) -> str:
    """Orient a reference window exactly as the Stage 02 matched background."""
    if strand == "sense":
        return reference_window
    if strand == "antisense":
        return reverse_complement(reference_window)
    raise ValueError(f"invalid expected strand: {strand}")


def physical_position(sequence: str, position_5p: int) -> str:
    if position_5p < 1 or position_5p > len(sequence):
        raise ValueError("physical position lies outside sequence")
    return sequence[position_5p - 1]


def position_from_3p(length: int, position_5p: int) -> int:
    return length - position_5p + 1


def gc9_14_fraction(sequence: str) -> float:
    if len(sequence) < 14:
        raise ValueError("GC9-14 requires a sequence of length at least 14")
    return sum(base in {"G", "C"} for base in sequence[8:14]) / 6.0


def safe_fraction(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator != 0 else None


def safe_delta(left: float | None, right: float | None) -> float | None:
    return None if left is None or right is None else left - right


def representation_enrichment(
    observed_fraction: float | None,
    expected_fraction: float | None,
    observed_total: float,
) -> float | None:
    if observed_total == 0 or observed_fraction is None or expected_fraction in (None, 0):
        return None
    return observed_fraction / expected_fraction


def accumulation_metrics(
    unique_fraction: float | None, abundance_fraction: float | None
) -> tuple[float | None, float | None, float | None]:
    delta = safe_delta(abundance_fraction, unique_fraction)
    if unique_fraction is None or abundance_fraction is None or unique_fraction <= 0:
        return None, delta, None
    ratio = abundance_fraction / unique_fraction
    log_ratio = math.log2(ratio) if math.isfinite(ratio) and ratio > 0 else None
    return ratio, delta, log_ratio


def finite(values: Iterable[float | None]) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(float(value))]


def median_or_none(values: Iterable[float | None]) -> float | None:
    usable = finite(values)
    return float(statistics.median(usable)) if usable else None


def percentile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return float(ordered[lower] * (1 - weight) + ordered[upper] * weight)


def sample_clustered_bootstrap(
    sample_values: dict[str, float],
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
    level: float = CI_LEVEL,
) -> tuple[float | None, float | None, int]:
    samples = sorted(sample_values)
    if not samples:
        return None, None, 0
    rng = random.Random(seed)
    estimates: list[float] = []
    for _ in range(replicates):
        selected = [rng.choice(samples) for _ in samples]
        estimates.append(float(statistics.median(sample_values[sample] for sample in selected)))
    alpha = (1 - level) / 2
    return percentile(estimates, alpha), percentile(estimates, 1 - alpha), len(estimates)


def exact_sign_test(values: Iterable[float | None]) -> dict[str, object]:
    usable = finite(values)
    positive = sum(value > 0 for value in usable)
    negative = sum(value < 0 for value in usable)
    nonzero = positive + negative
    if nonzero == 0:
        raw_p = None
        estimability = "not_estimable_no_nonzero_sample_deltas"
    else:
        tail = min(positive, negative)
        probability = sum(math.comb(nonzero, count) for count in range(tail + 1)) / (2 ** nonzero)
        raw_p = min(1.0, 2.0 * probability)
        estimability = "estimable"
    return {
        "n_total": len(usable),
        "n_nonzero": nonzero,
        "n_positive": positive,
        "n_negative": negative,
        "raw_p": raw_p,
        "estimability": estimability,
    }


def adjust_pvalues(values: Sequence[float | None], method: str) -> list[float | None]:
    indexed = [(index, float(value)) for index, value in enumerate(values) if value is not None]
    output: list[float | None] = [None] * len(values)
    if not indexed:
        return output
    indexed.sort(key=lambda item: item[1])
    count = len(indexed)
    factor = 1.0 if method == "BH" else sum(1.0 / rank for rank in range(1, count + 1))
    if method not in {"BH", "BY"}:
        raise ValueError(f"unsupported correction method: {method}")
    running = 1.0
    for rank in range(count, 0, -1):
        index, p_value = indexed[rank - 1]
        adjusted = min(1.0, p_value * count * factor / rank)
        running = min(running, adjusted)
        output[index] = running
    return output


def parse_fasta(path: Path) -> list[str]:
    records: list[str] = []
    sequence: list[str] = []
    seen_header = False
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    records.append("".join(sequence).upper())
                seen_header = True
                sequence = []
            else:
                if not seen_header:
                    raise Stage07Error(f"sequence before FASTA header in {path}:{line_number}")
                sequence.append(line)
    if seen_header:
        records.append("".join(sequence).upper())
    if not records or any(not record for record in records):
        raise Stage07Error(f"empty FASTA record in {path}")
    return records


def enumerate_background(records: Sequence[str], length: int) -> dict[str, object]:
    counts: dict[str, Counter] = {strand: Counter() for strand in STRANDS}
    gc_sums = {strand: 0.0 for strand in STRANDS}
    valid = candidate = 0
    for record in records:
        starts = max(0, len(record) - length + 1)
        candidate += starts
        for start in range(starts):
            window = record[start : start + length]
            if any(base not in NUCLEOTIDES for base in window):
                continue
            valid += 1
            oriented = {
                strand: expected_physical_sequence(window, strand) for strand in STRANDS
            }
            for strand, sequence in oriented.items():
                gc_sums[strand] += gc9_14_fraction(sequence)
                for position, base in enumerate(sequence, start=1):
                    counts[strand][(position, base)] += 1
    return {
        "counts": counts,
        "gc_sums": gc_sums,
        "valid": valid,
        "candidate": candidate,
        "excluded": candidate - valid,
    }


def iter_feature_rows(root: Path, samples: Iterable[str]) -> Iterator[dict[str, str]]:
    for sample in sorted(samples):
        path = root / "tables" / sample / f"{sample}.read_level_features.tsv.gz"
        with gzip.open(path, "rt", newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def _representation_endpoint(mode: str) -> str:
    return "unique_representation" if mode == "unique_sequence" else "abundance_representation"


def build_pair_tables(
    eligibility: Sequence[dict[str, str]],
    feature_rows: Iterable[dict[str, str]],
    backgrounds: dict[tuple[str, str], Sequence[str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    metadata = {
        (row["sample"], row["analysis_unit"]): row
        for row in eligibility
        if is_true(row.get("primary_eligible", ""))
    }
    pairs = set(metadata)
    primary_samples = {sample for sample, _unit in pairs}
    if len(primary_samples) != EXPECTED_PRIMARY_SAMPLES or len(pairs) != EXPECTED_PRIMARY_PAIRS:
        raise Stage07Error(
            "primary eligibility scope mismatch: "
            f"expected {EXPECTED_PRIMARY_SAMPLES} samples/{EXPECTED_PRIMARY_PAIRS} pairs, "
            f"observed {len(primary_samples)} samples/{len(pairs)} pairs"
        )
    if set(backgrounds) != pairs:
        raise Stage07Error("background sample-virus scope does not match primary eligibility")

    abundance_counts = {pair: Counter() for pair in pairs}
    abundance_totals = {pair: Counter() for pair in pairs}
    abundance_gc = {pair: Counter() for pair in pairs}
    unique_counts = {pair: Counter() for pair in pairs}
    unique_totals = {pair: Counter() for pair in pairs}
    unique_gc = {pair: Counter() for pair in pairs}
    sample_unique: dict[tuple[str, str, int, str], set[str]] = defaultdict(set)
    current_sample: str | None = None
    completed_samples: set[str] = set()
    rows_examined = retained_rows = length_mismatches = unexpected_bases = 0
    category_values: dict[str, set[str]] = defaultdict(set)

    def add_sequence(
        counts: Counter,
        totals: Counter,
        gc_totals: Counter,
        pair: tuple[str, str],
        length: int,
        strand: str,
        sequence: str,
        weight: float,
    ) -> None:
        totals[pair][(length, strand)] += weight
        gc_totals[pair][(length, strand)] += weight * gc9_14_fraction(sequence)
        for position, base in enumerate(sequence, start=1):
            counts[pair][(length, strand, position, base)] += weight

    def flush_unique() -> None:
        for (sample, unit, length, strand), sequences in sample_unique.items():
            pair = (sample, unit)
            for sequence in sequences:
                add_sequence(
                    unique_counts, unique_totals, unique_gc,
                    pair, length, strand, sequence, 1.0,
                )
        sample_unique.clear()

    for row in feature_rows:
        rows_examined += 1
        sample = row.get("sample", "")
        if current_sample is None:
            current_sample = sample
        elif sample != current_sample:
            flush_unique()
            completed_samples.add(current_sample)
            if sample in completed_samples:
                raise Stage07Error("feature rows must be grouped by sample")
            current_sample = sample
        for column in EXPECTED_CATEGORIES:
            category_values[column].add(row.get(column, ""))
        pair = (sample, row.get("virus", ""))
        if (
            pair not in pairs
            or row.get("mapping_mode") != "exact"
            or row.get("virus_assignment") != "assigned"
            or row.get("strand") not in STRANDS
        ):
            continue
        try:
            length = int(row["length"])
        except (KeyError, TypeError, ValueError) as exc:
            raise Stage07Error(f"invalid declared length at feature row {rows_examined}") from exc
        if length not in TARGET_LENGTHS:
            continue
        sequence = row.get("sequence", "")
        if len(sequence) != length:
            length_mismatches += 1
            continue
        if any(base not in NUCLEOTIDES for base in sequence):
            unexpected_bases += 1
            continue
        try:
            count = float(row["count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise Stage07Error(f"invalid count at feature row {rows_examined}") from exc
        if count < 0 or not math.isfinite(count):
            raise Stage07Error(f"invalid count at feature row {rows_examined}")
        retained_rows += 1
        strand = row["strand"]
        sequence = observed_physical_sequence(sequence, strand)
        add_sequence(
            abundance_counts, abundance_totals, abundance_gc,
            pair, length, strand, sequence, count,
        )
        sample_unique[(pair[0], pair[1], length, strand)].add(sequence)
    flush_unique()

    background_stats = {
        (pair, length): enumerate_background(backgrounds[pair], length)
        for pair in pairs
        for length in TARGET_LENGTHS
    }
    pair_rows: list[dict[str, object]] = []
    gc_rows: list[dict[str, object]] = []
    observed_deviations: list[float] = []
    expected_deviations: list[float] = []

    for pair in sorted(pairs):
        meta = metadata[pair]
        for length in TARGET_LENGTHS:
            stats = background_stats[(pair, length)]
            valid = int(stats["valid"])
            for strand in STRANDS:
                observed_lookup: dict[tuple[str, int, str], float | None] = {}
                for mode, counters, totals in (
                    ("unique_sequence", unique_counts, unique_totals),
                    ("abundance", abundance_counts, abundance_totals),
                ):
                    total = float(totals[pair][(length, strand)])
                    for position in range(1, length + 1):
                        fractions: list[float] = []
                        for nucleotide in NUCLEOTIDES:
                            weight = float(counters[pair][(length, strand, position, nucleotide)])
                            fraction = safe_fraction(weight, total)
                            observed_lookup[(mode, position, nucleotide)] = fraction
                            if fraction is not None:
                                fractions.append(fraction)
                        if fractions:
                            observed_deviations.append(abs(sum(fractions) - 1.0))

                for position in range(1, length + 1):
                    expected_at_position: list[float] = []
                    for nucleotide in NUCLEOTIDES:
                        expected_weight = float(stats["counts"][strand][(position, nucleotide)])
                        expected_fraction = safe_fraction(expected_weight, valid)
                        if expected_fraction is not None:
                            expected_at_position.append(expected_fraction)
                        unique_fraction = observed_lookup[("unique_sequence", position, nucleotide)]
                        abundance_fraction = observed_lookup[("abundance", position, nucleotide)]
                        accumulation_ratio, accumulation_delta, log_ratio = accumulation_metrics(
                            unique_fraction, abundance_fraction
                        )
                        common = {
                            "sample": pair[0],
                            "analysis_unit": pair[1],
                            "biological_virus": meta["biological_virus"],
                            "polarity": meta["polarity"],
                            "length": length,
                            "strand": strand,
                            "position_5p": position,
                            "position_from_3p": position_from_3p(length, position),
                            "nucleotide": nucleotide,
                            "expected_nucleotide_weight": expected_weight if valid else None,
                            "expected_total_windows": valid if valid else None,
                            "expected_fraction": expected_fraction,
                            "unique_fraction": unique_fraction,
                            "abundance_fraction": abundance_fraction,
                            "accumulation_ratio": accumulation_ratio,
                            "accumulation_delta_fraction": accumulation_delta,
                            "log2_accumulation_ratio": log_ratio,
                            "valid_background_windows": valid,
                        }
                        for mode, counters, totals in (
                            ("unique_sequence", unique_counts, unique_totals),
                            ("abundance", abundance_counts, abundance_totals),
                        ):
                            total = float(totals[pair][(length, strand)])
                            observed_weight = float(
                                counters[pair][(length, strand, position, nucleotide)]
                            )
                            observed_fraction = observed_lookup[(mode, position, nucleotide)]
                            pair_rows.append(
                                common
                                | {
                                    "weighting_mode": mode,
                                    "observed_nucleotide_weight": observed_weight,
                                    "observed_total_weight": total,
                                    "observed_fraction": observed_fraction,
                                    "representation_enrichment": representation_enrichment(
                                        observed_fraction, expected_fraction, total
                                    ),
                                    "representation_delta_fraction": safe_delta(
                                        observed_fraction, expected_fraction
                                    ),
                                }
                            )
                        pair_rows.append(
                            common
                            | {
                                "weighting_mode": "accumulation",
                                "observed_nucleotide_weight": None,
                                "observed_total_weight": None,
                                "observed_fraction": None,
                                "representation_enrichment": None,
                                "representation_delta_fraction": None,
                            }
                        )
                    if expected_at_position:
                        expected_deviations.append(abs(sum(expected_at_position) - 1.0))

                unique_total = float(unique_totals[pair][(length, strand)])
                abundance_total = float(abundance_totals[pair][(length, strand)])
                expected_gc = safe_fraction(float(stats["gc_sums"][strand]), valid)
                unique_gc_mean = safe_fraction(
                    float(unique_gc[pair][(length, strand)]), unique_total
                )
                abundance_gc_mean = safe_fraction(
                    float(abundance_gc[pair][(length, strand)]), abundance_total
                )
                gc_rows.append(
                    {
                        "sample": pair[0],
                        "analysis_unit": pair[1],
                        "biological_virus": meta["biological_virus"],
                        "polarity": meta["polarity"],
                        "length": length,
                        "strand": strand,
                        "observed_GC9_14_mean_unique": unique_gc_mean,
                        "observed_GC9_14_mean_abundance": abundance_gc_mean,
                        "expected_GC9_14_mean": expected_gc,
                        "GC9_14_delta_unique_vs_expected": safe_delta(unique_gc_mean, expected_gc),
                        "GC9_14_delta_abundance_vs_expected": safe_delta(
                            abundance_gc_mean, expected_gc
                        ),
                        "GC9_14_accumulation_delta": safe_delta(
                            abundance_gc_mean, unique_gc_mean
                        ),
                    }
                )

    stats_out: dict[str, object] = {
        "primary_samples": len(primary_samples),
        "primary_pairs": len(pairs),
        "rows_examined": rows_examined,
        "retained_rows": retained_rows,
        "length_mismatches": length_mismatches,
        "unexpected_bases": unexpected_bases,
        "category_values": category_values,
        "background_stats": background_stats,
        "background_records": sum(len(records) for records in backgrounds.values()),
        "max_observed_deviation": max(observed_deviations, default=0.0),
        "max_expected_deviation": max(expected_deviations, default=0.0),
    }
    return pair_rows, gc_rows, stats_out


def aggregate_positional(
    pair_rows: Sequence[dict[str, object]],
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        key = tuple(
            row[field]
            for field in ("length", "strand", "weighting_mode", "position_5p", "nucleotide")
        )
        grouped[key].append(row)
    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    metrics = {
        "observed_fraction": "sample_observed_fraction",
        "expected_fraction": "sample_expected_fraction",
        "representation_enrichment": "sample_representation_enrichment",
        "representation_delta_fraction": "sample_representation_delta_fraction",
        "unique_fraction": "sample_unique_fraction",
        "abundance_fraction": "sample_abundance_fraction",
        "accumulation_ratio": "sample_accumulation_ratio",
        "accumulation_delta_fraction": "sample_accumulation_delta_fraction",
        "log2_accumulation_ratio": "sample_log2_accumulation_ratio",
    }

    for key in sorted(grouped):
        length, strand, mode, position, nucleotide = key
        rows = grouped[key]
        by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            by_sample[str(row["sample"])].append(row)
        current_samples: list[dict[str, object]] = []
        for sample in sorted(by_sample):
            values = by_sample[sample]
            sample_row: dict[str, object] = {
                "sample": sample,
                "length": length,
                "strand": strand,
                "weighting_mode": mode,
                "position_5p": position,
                "position_from_3p": position_from_3p(int(length), int(position)),
                "nucleotide": nucleotide,
            }
            for pair_field, sample_field in metrics.items():
                sample_row[sample_field] = median_or_none(row[pair_field] for row in values)
            inference_field = (
                "sample_accumulation_delta_fraction"
                if mode == "accumulation"
                else "sample_representation_delta_fraction"
            )
            sample_row["n_virus_units"] = sum(
                value is not None and math.isfinite(float(value))
                for value in (row[
                    "accumulation_delta_fraction"
                    if mode == "accumulation"
                    else "representation_delta_fraction"
                ] for row in values)
            )
            if sample_row[inference_field] is not None:
                current_samples.append(sample_row)
            sample_rows.append(sample_row)

        summary: dict[str, object] = {
            "strand": strand,
            "length": length,
            "position_5p": position,
            "position_from_3p": position_from_3p(int(length), int(position)),
            "nucleotide": nucleotide,
            "weighting_mode": mode,
            "endpoint": "accumulation" if mode == "accumulation" else _representation_endpoint(str(mode)),
            "n_samples_total": len(current_samples),
        }
        summary_mapping = {
            "sample_observed_fraction": "sample_balanced_observed_fraction",
            "sample_expected_fraction": "sample_balanced_expected_fraction",
            "sample_representation_enrichment": "sample_balanced_representation_enrichment",
            "sample_representation_delta_fraction": "sample_balanced_representation_delta_fraction",
            "sample_unique_fraction": "sample_balanced_unique_fraction",
            "sample_abundance_fraction": "sample_balanced_abundance_fraction",
            "sample_accumulation_ratio": "sample_balanced_accumulation_ratio",
            "sample_accumulation_delta_fraction": "sample_balanced_accumulation_delta_fraction",
            "sample_log2_accumulation_ratio": "log2_accumulation_ratio",
        }
        for sample_field, summary_field in summary_mapping.items():
            summary[summary_field] = median_or_none(row[sample_field] for row in current_samples)

        effect_field = (
            "sample_accumulation_ratio"
            if mode == "accumulation"
            else "sample_representation_enrichment"
        )
        effect_values = {
            str(row["sample"]): float(row[effect_field])
            for row in current_samples
            if row[effect_field] is not None and math.isfinite(float(row[effect_field]))
        }
        ci_low, ci_high, valid = sample_clustered_bootstrap(
            effect_values, bootstrap_replicates, seed, CI_LEVEL
        )
        delta_field = (
            "sample_accumulation_delta_fraction"
            if mode == "accumulation"
            else "sample_representation_delta_fraction"
        )
        sign = exact_sign_test(row[delta_field] for row in current_samples)
        summary.update(
            {
                "n_samples": len(effect_values),
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "bootstrap_replicates_requested": bootstrap_replicates,
                "bootstrap_replicates_valid": valid,
                "bootstrap_seed": seed,
                "ci_method": "percentile",
                "ci_level": CI_LEVEL,
                "sign_test_n_nonzero": sign["n_nonzero"],
                "sign_test_n_positive": sign["n_positive"],
                "sign_test_n_negative": sign["n_negative"],
                "sign_test_estimability": sign["estimability"],
                "raw_p": sign["raw_p"],
                "bh_p": None,
                "by_p": None,
            }
        )
        summary_rows.append(summary)

    apply_discovery_corrections(summary_rows)
    return sample_rows, summary_rows


def apply_discovery_corrections(rows: list[dict[str, object]]) -> None:
    families: dict[tuple[object, ...], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        length = int(row["length"])
        position = int(row["position_5p"])
        if 3 <= position <= length - 2:
            families[(length, row["endpoint"], row["strand"])].append(index)
    for indices in families.values():
        raw = [rows[index]["raw_p"] for index in indices]
        bh = adjust_pvalues(raw, "BH")
        by = adjust_pvalues(raw, "BY")
        for index, bh_value, by_value in zip(indices, bh, by):
            rows[index]["bh_p"] = bh_value
            rows[index]["by_p"] = by_value


def aggregate_gc(
    pair_rows: Sequence[dict[str, object]],
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        grouped[(row["length"], row["strand"])].append(row)
    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    delta_fields = {
        "unique_representation": "GC9_14_delta_unique_vs_expected",
        "abundance_representation": "GC9_14_delta_abundance_vs_expected",
        "accumulation": "GC9_14_accumulation_delta",
    }
    for (length, strand), rows in sorted(grouped.items()):
        by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            by_sample[str(row["sample"])].append(row)
        current_samples: list[dict[str, object]] = []
        for sample in sorted(by_sample):
            values = by_sample[sample]
            sample_row: dict[str, object] = {
                "sample": sample,
                "length": length,
                "strand": strand,
                "n_virus_units": len(values),
            }
            for field in GC_PAIR_FIELDS[6:]:
                sample_row[field] = median_or_none(row[field] for row in values)
            current_samples.append(sample_row)
            sample_rows.append(sample_row)
        for endpoint, effect_field in delta_fields.items():
            sample_values = {
                str(row["sample"]): float(row[effect_field])
                for row in current_samples
                if row[effect_field] is not None and math.isfinite(float(row[effect_field]))
            }
            estimate = median_or_none(sample_values.values())
            ci_low, ci_high, valid = sample_clustered_bootstrap(
                sample_values, bootstrap_replicates, seed, CI_LEVEL
            )
            sign = exact_sign_test(sample_values.values())
            summary_rows.append(
                {
                    "strand": strand,
                    "length": length,
                    "endpoint": endpoint,
                    "effect_metric": effect_field,
                    "n_samples": len(sample_values),
                    "effect_estimate": estimate,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "bootstrap_replicates_requested": bootstrap_replicates,
                    "bootstrap_replicates_valid": valid,
                    "bootstrap_seed": seed,
                    "ci_method": "percentile",
                    "ci_level": CI_LEVEL,
                    "sign_test_n_nonzero": sign["n_nonzero"],
                    "sign_test_n_positive": sign["n_positive"],
                    "sign_test_n_negative": sign["n_negative"],
                    "sign_test_estimability": sign["estimability"],
                    "raw_p": sign["raw_p"],
                    "bh_p": None,
                }
            )
    return sample_rows, summary_rows


def build_literature_validation(
    positional_summary: Sequence[dict[str, object]],
    gc_summary: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in positional_summary:
        if (
            row["strand"] == "antisense"
            and int(row["position_5p"]) == 10
            and row["nucleotide"] == "A"
        ):
            mode = str(row["weighting_mode"])
            if mode == "accumulation":
                effect_metric = "sample_balanced_accumulation_ratio"
                effect = row[effect_metric]
                delta = row["sample_balanced_accumulation_delta_fraction"]
            else:
                effect_metric = "sample_balanced_representation_enrichment"
                effect = row[effect_metric]
                delta = row["sample_balanced_representation_delta_fraction"]
            output.append(
                {
                    "feature": "A10",
                    "strand": "antisense",
                    "length": row["length"],
                    "endpoint": row["endpoint"],
                    "effect_metric": effect_metric,
                    "effect_estimate": effect,
                    "paired_delta_estimate": delta,
                    "bootstrap_ci_low": row["bootstrap_ci_low"],
                    "bootstrap_ci_high": row["bootstrap_ci_high"],
                    "n_samples": row["n_samples_total"],
                    "sign_test_n_nonzero": row["sign_test_n_nonzero"],
                    "sign_test_n_positive": row["sign_test_n_positive"],
                    "sign_test_n_negative": row["sign_test_n_negative"],
                    "sign_test_estimability": row["sign_test_estimability"],
                    "raw_p": row["raw_p"],
                    "bh_p": None,
                }
            )
    for row in gc_summary:
        if row["strand"] == "antisense":
            output.append(
                {
                    "feature": "GC9_14_continuous",
                    "strand": "antisense",
                    "length": row["length"],
                    "endpoint": row["endpoint"],
                    "effect_metric": row["effect_metric"],
                    "effect_estimate": row["effect_estimate"],
                    "paired_delta_estimate": row["effect_estimate"],
                    "bootstrap_ci_low": row["bootstrap_ci_low"],
                    "bootstrap_ci_high": row["bootstrap_ci_high"],
                    "n_samples": row["n_samples"],
                    "sign_test_n_nonzero": row["sign_test_n_nonzero"],
                    "sign_test_n_positive": row["sign_test_n_positive"],
                    "sign_test_n_negative": row["sign_test_n_negative"],
                    "sign_test_estimability": row["sign_test_estimability"],
                    "raw_p": row["raw_p"],
                    "bh_p": None,
                }
            )
    output.sort(key=lambda row: (str(row["feature"]), int(row["length"]), str(row["endpoint"])))
    adjusted = adjust_pvalues([row["raw_p"] for row in output], "BH")
    for row, value in zip(output, adjusted):
        row["bh_p"] = value
    return output


def regional_windows(length: int, width: int = REGIONAL_GC_WIDTH) -> list[dict[str, object]]:
    """Enumerate the canonical fixed-width regional-GC coordinates."""
    if width != REGIONAL_GC_WIDTH:
        raise ValueError(f"canonical regional GC width must be {REGIONAL_GC_WIDTH}")
    if length < width:
        raise ValueError("regional GC width exceeds sequence length")
    output: list[dict[str, object]] = []
    for start in range(1, length - width + 2):
        end = start + width - 1
        near = length - end + 1
        far = length - start + 1
        output.append(
            {
                "start_5p": start,
                "end_5p": end,
                "near_3p": near,
                "far_3p": far,
                "region_5p": f"GC{start}-{end}",
                "region_3p": f"GC_3p{near}-{far}",
            }
        )
    return output


def regional_gc6_fraction(sequence: str, start_5p: int) -> float:
    """Direct six-position GC fraction for deterministic regression tests."""
    end_5p = start_5p + REGIONAL_GC_WIDTH - 1
    if start_5p < 1 or end_5p > len(sequence):
        raise ValueError("regional GC window lies outside sequence")
    region = sequence[start_5p - 1 : end_5p]
    return sum(base in {"G", "C"} for base in region) / REGIONAL_GC_WIDTH


def regional_gc_from_positional_fractions(
    gc_fraction_by_position: dict[int, float | None], start_5p: int
) -> float | None:
    """Derive mean regional GC from constituent positional P(G)+P(C)."""
    values = [
        gc_fraction_by_position.get(position)
        for position in range(start_5p, start_5p + REGIONAL_GC_WIDTH)
    ]
    if any(value is None for value in values):
        return None
    return sum(float(value) for value in values) / REGIONAL_GC_WIDTH


def regional_gc_deltas(
    observed_unique: float | None,
    observed_abundance: float | None,
    expected: float | None,
) -> tuple[float | None, float | None, float | None]:
    return (
        safe_delta(observed_unique, expected),
        safe_delta(observed_abundance, expected),
        safe_delta(observed_abundance, observed_unique),
    )


def build_regional_gc_pair_rows(
    positional_pair_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    """Derive all six-nt regional means from existing positional fractions."""
    representation_rows = [
        row for row in positional_pair_rows if row["weighting_mode"] in REPRESENTATION_MODES
    ]
    identities: dict[tuple[object, ...], dict[str, object]] = {}
    lookup: dict[tuple[object, ...], dict[str, object]] = {}
    for row in representation_rows:
        identity = (
            row["sample"], row["analysis_unit"], int(row["length"]), row["strand"]
        )
        identities.setdefault(identity, row)
        lookup[
            identity
            + (row["weighting_mode"], int(row["position_5p"]), row["nucleotide"])
        ] = row

    output: list[dict[str, object]] = []
    for identity in sorted(identities):
        sample, unit, length, strand = identity
        representative = identities[identity]
        positional_gc: dict[str, dict[int, float | None]] = {
            "unique_sequence": {}, "abundance": {}, "expected": {}
        }
        for position in range(1, int(length) + 1):
            for mode in REPRESENTATION_MODES:
                g_row = lookup[identity + (mode, position, "G")]
                c_row = lookup[identity + (mode, position, "C")]
                g_fraction = g_row["observed_fraction"]
                c_fraction = c_row["observed_fraction"]
                positional_gc[mode][position] = (
                    None
                    if g_fraction is None or c_fraction is None
                    else float(g_fraction) + float(c_fraction)
                )
            expected_g = lookup[identity + ("unique_sequence", position, "G")][
                "expected_fraction"
            ]
            expected_c = lookup[identity + ("unique_sequence", position, "C")][
                "expected_fraction"
            ]
            positional_gc["expected"][position] = (
                None
                if expected_g is None or expected_c is None
                else float(expected_g) + float(expected_c)
            )

        for coordinates in regional_windows(int(length)):
            start = int(coordinates["start_5p"])
            unique_mean = regional_gc_from_positional_fractions(
                positional_gc["unique_sequence"], start
            )
            abundance_mean = regional_gc_from_positional_fractions(
                positional_gc["abundance"], start
            )
            expected_mean = regional_gc_from_positional_fractions(
                positional_gc["expected"], start
            )
            unique_delta, abundance_delta, accumulation_delta = regional_gc_deltas(
                unique_mean, abundance_mean, expected_mean
            )
            output.append(
                {
                    "sample": sample,
                    "analysis_unit": unit,
                    "biological_virus": representative["biological_virus"],
                    "polarity": representative["polarity"],
                    "strand": strand,
                    "length": length,
                    **coordinates,
                    "observed_gc6_mean_unique": unique_mean,
                    "observed_gc6_mean_abundance": abundance_mean,
                    "expected_gc6_mean": expected_mean,
                    "regional_gc6_delta_unique_vs_expected": unique_delta,
                    "regional_gc6_delta_abundance_vs_expected": abundance_delta,
                    "regional_gc6_accumulation_delta": accumulation_delta,
                }
            )
    return output


def apply_regional_gc_corrections(rows: list[dict[str, object]]) -> None:
    """Correct exploratory windows separately by length, endpoint, and strand."""
    families: dict[tuple[object, ...], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if row["evidence_class"] == "exploratory_regional_gc6":
            families[(row["length"], row["endpoint"], row["strand"])].append(index)
    for indices in families.values():
        raw = [rows[index]["raw_p"] for index in indices]
        bh = adjust_pvalues(raw, "BH")
        by = adjust_pvalues(raw, "BY")
        for index, bh_value, by_value in zip(indices, bh, by):
            rows[index]["regional_bh_p"] = bh_value
            rows[index]["regional_by_p"] = by_value


def aggregate_regional_gc(
    pair_rows: Sequence[dict[str, object]],
    literature_validation: Sequence[dict[str, object]],
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Reuse the existing GC aggregation/inference path for every six-nt window."""
    grouped: dict[tuple[int, int], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        grouped[(int(row["length"]), int(row["start_5p"]))].append(row)

    validation_bh = {
        (int(row["length"]), str(row["endpoint"])): row["bh_p"]
        for row in literature_validation
        if row["feature"] == "GC9_14_continuous" and row["strand"] == "antisense"
    }
    sample_output: list[dict[str, object]] = []
    summary_output: list[dict[str, object]] = []
    for (length, start), rows in sorted(grouped.items()):
        coordinates = regional_windows(length)[start - 1]
        transformed = [
            {
                "sample": row["sample"],
                "analysis_unit": row["analysis_unit"],
                "biological_virus": row["biological_virus"],
                "polarity": row["polarity"],
                "length": row["length"],
                "strand": row["strand"],
                "observed_GC9_14_mean_unique": row["observed_gc6_mean_unique"],
                "observed_GC9_14_mean_abundance": row["observed_gc6_mean_abundance"],
                "expected_GC9_14_mean": row["expected_gc6_mean"],
                "GC9_14_delta_unique_vs_expected": row[
                    "regional_gc6_delta_unique_vs_expected"
                ],
                "GC9_14_delta_abundance_vs_expected": row[
                    "regional_gc6_delta_abundance_vs_expected"
                ],
                "GC9_14_accumulation_delta": row["regional_gc6_accumulation_delta"],
            }
            for row in rows
        ]
        gc_samples, gc_summaries = aggregate_gc(transformed, bootstrap_replicates, seed)
        regional_samples: list[dict[str, object]] = []
        for row in gc_samples:
            regional_row = {
                "sample": row["sample"],
                "strand": row["strand"],
                "length": row["length"],
                **coordinates,
                "n_virus_units": row["n_virus_units"],
                "observed_gc6_mean_unique": row["observed_GC9_14_mean_unique"],
                "observed_gc6_mean_abundance": row["observed_GC9_14_mean_abundance"],
                "expected_gc6_mean": row["expected_GC9_14_mean"],
                "regional_gc6_delta_unique_vs_expected": row[
                    "GC9_14_delta_unique_vs_expected"
                ],
                "regional_gc6_delta_abundance_vs_expected": row[
                    "GC9_14_delta_abundance_vs_expected"
                ],
                "regional_gc6_accumulation_delta": row["GC9_14_accumulation_delta"],
            }
            regional_samples.append(regional_row)
            sample_output.append(regional_row)

        for row in gc_summaries:
            endpoint = str(row["endpoint"])
            strand = str(row["strand"])
            evidence_class = (
                "literature_validation_gc9_14"
                if start == 9
                else "exploratory_regional_gc6"
            )
            relevant_samples = [item for item in regional_samples if item["strand"] == strand]
            if endpoint == "unique_representation":
                observed_field = "observed_gc6_mean_unique"
                expected_field = "expected_gc6_mean"
            elif endpoint == "abundance_representation":
                observed_field = "observed_gc6_mean_abundance"
                expected_field = "expected_gc6_mean"
            else:
                observed_field = expected_field = None
            summary_output.append(
                {
                    "strand": strand,
                    "length": length,
                    **coordinates,
                    "endpoint": endpoint,
                    "n_samples": row["n_samples"],
                    "sample_balanced_observed_gc6": (
                        None
                        if observed_field is None
                        else median_or_none(item[observed_field] for item in relevant_samples)
                    ),
                    "sample_balanced_expected_gc6": (
                        None
                        if expected_field is None
                        else median_or_none(item[expected_field] for item in relevant_samples)
                    ),
                    "sample_balanced_regional_gc6_delta": row["effect_estimate"],
                    "bootstrap_ci_low": row["bootstrap_ci_low"],
                    "bootstrap_ci_high": row["bootstrap_ci_high"],
                    "bootstrap_replicates_requested": row["bootstrap_replicates_requested"],
                    "bootstrap_replicates_valid": row["bootstrap_replicates_valid"],
                    "bootstrap_seed": row["bootstrap_seed"],
                    "ci_method": row["ci_method"],
                    "ci_level": row["ci_level"],
                    "sign_test_n_nonzero": row["sign_test_n_nonzero"],
                    "sign_test_n_positive": row["sign_test_n_positive"],
                    "sign_test_n_negative": row["sign_test_n_negative"],
                    "sign_test_estimability": row["sign_test_estimability"],
                    "raw_p": row["raw_p"],
                    "evidence_class": evidence_class,
                    "validation_bh_p": (
                        validation_bh.get((length, endpoint))
                        if start == 9 and strand == "antisense"
                        else None
                    ),
                    "regional_bh_p": None,
                    "regional_by_p": None,
                }
            )

    apply_regional_gc_corrections(summary_output)

    summary_output.sort(
        key=lambda row: (
            int(row["length"]), str(row["strand"]), int(row["start_5p"]),
            str(row["endpoint"]),
        )
    )
    sample_output.sort(
        key=lambda row: (
            str(row["sample"]), int(row["length"]), str(row["strand"]),
            int(row["start_5p"]),
        )
    )
    discovery = [
        row for row in summary_output if row["evidence_class"] == "exploratory_regional_gc6"
    ]
    return sample_output, summary_output, discovery


def regional_gc9_14_regression(
    regional_pair: Sequence[dict[str, object]],
    regional_sample: Sequence[dict[str, object]],
    regional_summary: Sequence[dict[str, object]],
    gc_pair: Sequence[dict[str, object]],
    gc_sample: Sequence[dict[str, object]],
    gc_summary: Sequence[dict[str, object]],
) -> tuple[bool, float, int]:
    """Require exact agreement of regional GC9-14 with the existing implementation."""
    differences: list[float] = []

    def compare(left: object, right: object) -> None:
        if left is None and right is None:
            differences.append(0.0)
        elif left is None or right is None:
            differences.append(math.inf)
        else:
            differences.append(abs(float(left) - float(right)))

    regional_pair_lookup = {
        (row["sample"], row["analysis_unit"], row["length"], row["strand"]): row
        for row in regional_pair if int(row["start_5p"]) == 9
    }
    pair_fields = {
        "observed_gc6_mean_unique": "observed_GC9_14_mean_unique",
        "observed_gc6_mean_abundance": "observed_GC9_14_mean_abundance",
        "expected_gc6_mean": "expected_GC9_14_mean",
        "regional_gc6_delta_unique_vs_expected": "GC9_14_delta_unique_vs_expected",
        "regional_gc6_delta_abundance_vs_expected": "GC9_14_delta_abundance_vs_expected",
        "regional_gc6_accumulation_delta": "GC9_14_accumulation_delta",
    }
    for old in gc_pair:
        new = regional_pair_lookup[(old["sample"], old["analysis_unit"], old["length"], old["strand"])]
        for new_field, old_field in pair_fields.items():
            compare(new[new_field], old[old_field])

    regional_sample_lookup = {
        (row["sample"], row["length"], row["strand"]): row
        for row in regional_sample if int(row["start_5p"]) == 9
    }
    for old in gc_sample:
        new = regional_sample_lookup[(old["sample"], old["length"], old["strand"])]
        for new_field, old_field in pair_fields.items():
            compare(new[new_field], old[old_field])

    regional_summary_lookup = {
        (row["length"], row["strand"], row["endpoint"]): row
        for row in regional_summary if int(row["start_5p"]) == 9
    }
    for old in gc_summary:
        new = regional_summary_lookup[(old["length"], old["strand"], old["endpoint"])]
        compare(new["sample_balanced_regional_gc6_delta"], old["effect_estimate"])

    maximum = max(differences, default=0.0)
    return maximum <= REGRESSION_TOLERANCE, maximum, len(differences)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_optional_float(value: object) -> float | None:
    text = str(value).strip()
    if text in {"", "NA", "nan", "None"}:
        return None
    parsed = float(text)
    return parsed if math.isfinite(parsed) else None


def _regression_difference(stage02_value: float | None, stage07_value: float | None) -> float | None:
    if stage02_value is None and stage07_value is None:
        return 0.0
    if stage02_value is None or stage07_value is None:
        return None
    stage07_formatted = float(format(stage07_value, ".12g"))
    return abs(stage02_value - stage07_formatted)


def stage02_terminal_regression(
    positional_by_pair: Sequence[dict[str, object]],
    positional_summary: Sequence[dict[str, object]],
    stage02_pair_path: Path,
    stage02_across_path: Path,
) -> tuple[list[dict[str, object]], bool, float, int]:
    terminal_mapping = {"5p1": 1, "5p2": 2, "3p2": -1, "3p1": 0}
    pair_lookup = {
        (
            row["sample"], row["analysis_unit"], int(row["length"]), row["strand"],
            row["weighting_mode"], int(row["position_5p"]), row["nucleotide"],
        ): row
        for row in positional_by_pair
        if row["weighting_mode"] in REPRESENTATION_MODES
    }
    summary_lookup = {
        (
            int(row["length"]), row["strand"], row["weighting_mode"],
            int(row["position_5p"]), row["nucleotide"],
        ): row
        for row in positional_summary
        if row["weighting_mode"] in REPRESENTATION_MODES
    }
    checks: list[dict[str, object]] = []

    for stage02_row in read_tsv(stage02_pair_path):
        strand = stage02_row["strand_scope"]
        if strand not in STRANDS:
            continue
        length = int(stage02_row["length"])
        terminal = stage02_row["terminal_position"]
        mapped = terminal_mapping[terminal]
        position = length + mapped if mapped <= 0 else mapped
        key = (
            stage02_row["sample"], stage02_row["analysis_unit"], length, strand,
            stage02_row["weighting_mode"], position, stage02_row["nucleotide"],
        )
        stage07_row = pair_lookup.get(key)
        if stage07_row is None:
            raise Stage07Error(f"Stage 07 terminal pair row missing for key {key}")
        for stage02_field, stage07_field, metric in (
            ("observed_fraction", "observed_fraction", "observed_fraction"),
            ("expected_fraction", "expected_fraction", "expected_fraction"),
            ("enrichment_ratio", "representation_enrichment", "enrichment_ratio"),
        ):
            old = parse_optional_float(stage02_row[stage02_field])
            new = stage07_row[stage07_field]
            new_float = None if new is None else float(new)
            difference = _regression_difference(old, new_float)
            status = "PASS" if difference is not None and difference <= REGRESSION_TOLERANCE else "FAIL"
            checks.append(
                {
                    "level": "pair",
                    "sample": key[0],
                    "analysis_unit": key[1],
                    "length": length,
                    "strand": strand,
                    "weighting_mode": key[4],
                    "terminal_position": terminal,
                    "position_5p": position,
                    "nucleotide": key[6],
                    "metric": metric,
                    "stage02_value": old,
                    "stage07_value": new_float,
                    "absolute_difference": difference,
                    "status": status,
                }
            )

    for stage02_row in read_tsv(stage02_across_path):
        strand = stage02_row["strand_scope"]
        if strand not in STRANDS:
            continue
        length = int(stage02_row["length"])
        terminal = stage02_row["terminal_position"]
        mapped = terminal_mapping[terminal]
        position = length + mapped if mapped <= 0 else mapped
        key = (
            length, strand, stage02_row["weighting_mode"], position, stage02_row["nucleotide"]
        )
        stage07_row = summary_lookup.get(key)
        if stage07_row is None:
            raise Stage07Error(f"Stage 07 terminal summary row missing for key {key}")
        old = parse_optional_float(stage02_row["sample_balanced_median_enrichment_ratio"])
        new_value = stage07_row["sample_balanced_representation_enrichment"]
        new = None if new_value is None else float(new_value)
        difference = _regression_difference(old, new)
        status = "PASS" if difference is not None and difference <= REGRESSION_TOLERANCE else "FAIL"
        checks.append(
            {
                "level": "sample_balanced",
                "sample": "ALL",
                "analysis_unit": "ALL",
                "length": length,
                "strand": strand,
                "weighting_mode": key[2],
                "terminal_position": terminal,
                "position_5p": position,
                "nucleotide": key[4],
                "metric": "sample_balanced_enrichment_ratio",
                "stage02_value": old,
                "stage07_value": new,
                "absolute_difference": difference,
                "status": status,
            }
        )

    failed = any(row["status"] == "FAIL" for row in checks)
    differences = [float(row["absolute_difference"]) for row in checks if row["absolute_difference"] is not None]
    return checks, not failed, max(differences, default=0.0), len(checks)


def build_qc(
    stats: dict[str, object],
    regression_pass: bool,
    regression_max: float,
    regression_n: int,
    regional_summary: Sequence[dict[str, object]],
    regional_discovery: Sequence[dict[str, object]],
    regional_regression_pass: bool,
    regional_regression_max: float,
    regional_regression_n: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(metric: str, status: str, value: object, details: str = "") -> None:
        rows.append({"metric": metric, "status": status, "value": value, "details": details})

    add("primary_samples", "PASS" if stats["primary_samples"] == EXPECTED_PRIMARY_SAMPLES else "FAIL", stats["primary_samples"])
    add("primary_eligible_sample_virus_units", "PASS" if stats["primary_pairs"] == EXPECTED_PRIMARY_PAIRS else "FAIL", stats["primary_pairs"])
    add("read_level_rows_examined", "INFO", stats["rows_examined"])
    add("exact_assigned_eligible_23_24_rows_retained", "INFO", stats["retained_rows"])
    add("declared_length_sequence_length_mismatches", "FAIL" if stats["length_mismatches"] else "PASS", stats["length_mismatches"])
    add("unexpected_observed_bases", "FAIL" if stats["unexpected_bases"] else "PASS", stats["unexpected_bases"])
    category_values = stats["category_values"]
    for column, expected in EXPECTED_CATEGORIES.items():
        unexpected = sorted(category_values[column] - expected)
        add(
            f"unexpected_{column}",
            "WARN" if unexpected else "PASS",
            len(unexpected),
            f"observed={sorted(category_values[column])}; unexpected={unexpected}",
        )
    add("background_fasta_records_examined", "INFO", stats["background_records"])
    background_stats = stats["background_stats"]
    for length in TARGET_LENGTHS:
        add(
            f"valid_{length}nt_background_windows",
            "INFO",
            sum(int(value["valid"]) for (pair, item_length), value in background_stats.items() if item_length == length),
        )
        add(
            f"excluded_{length}nt_background_windows",
            "INFO",
            sum(int(value["excluded"]) for (pair, item_length), value in background_stats.items() if item_length == length),
        )
    zero_background = [key for key, value in background_stats.items() if int(value["valid"]) == 0]
    add("units_with_zero_valid_background_windows", "WARN" if zero_background else "PASS", len(zero_background), str(zero_background))
    add("maximum_observed_frequency_sum_deviation", "PASS" if stats["max_observed_deviation"] <= FREQUENCY_TOLERANCE else "FAIL", stats["max_observed_deviation"])
    add("maximum_expected_frequency_sum_deviation", "PASS" if stats["max_expected_deviation"] <= FREQUENCY_TOLERANCE else "FAIL", stats["max_expected_deviation"])
    add("all_physical_positions_enumerated", "PASS", sum(TARGET_LENGTHS), "all positions emitted for both lengths where populations are represented")
    add("pseudocounts_added", "PASS", 0)
    add("stage02_terminal_regression", "PASS" if regression_pass else "FAIL", regression_n, f"maximum_absolute_difference={regression_max:.17g}; tolerance={REGRESSION_TOLERANCE}")
    window_counts = {
        length: len({int(row["start_5p"]) for row in regional_summary if int(row["length"]) == length})
        for length in TARGET_LENGTHS
    }
    for length, expected in ((23, 18), (24, 19)):
        add(
            f"regional_gc6_{length}nt_window_count",
            "PASS" if window_counts[length] == expected else "FAIL",
            window_counts[length],
            f"expected={expected}",
        )
    widths = {
        int(row["end_5p"]) - int(row["start_5p"]) + 1 for row in regional_summary
    }
    add("regional_gc_window_widths", "PASS" if widths == {REGIONAL_GC_WIDTH} else "FAIL", ",".join(map(str, sorted(widths))))
    coordinate_ok = all(
        int(row["near_3p"]) == int(row["length"]) - int(row["end_5p"]) + 1
        and int(row["far_3p"]) == int(row["length"]) - int(row["start_5p"]) + 1
        for row in regional_summary
    )
    add("regional_gc_coordinate_conversion", "PASS" if coordinate_ok else "FAIL", int(coordinate_ok))
    coverage: dict[tuple[object, ...], set[str]] = defaultdict(set)
    for row in regional_summary:
        coverage[(row["length"], row["start_5p"], row["endpoint"])].add(str(row["strand"]))
    coverage_ok = all(strands == set(STRANDS) for strands in coverage.values())
    add("regional_gc_sense_antisense_coverage", "PASS" if coverage_ok else "FAIL", len(coverage))
    gc9_rows = [row for row in regional_summary if int(row["start_5p"]) == 9]
    gc9_excluded = (
        all(row["regional_bh_p"] is None and row["regional_by_p"] is None for row in gc9_rows)
        and all(int(row["start_5p"]) != 9 for row in regional_discovery)
    )
    add("regional_gc9_14_excluded_from_exploratory_families", "PASS" if gc9_excluded else "FAIL", int(gc9_excluded))
    family_counts: dict[tuple[object, ...], int] = Counter(
        (row["length"], row["endpoint"], row["strand"]) for row in regional_discovery
    )
    for length, expected in ((23, 17), (24, 18)):
        counts = [count for key, count in family_counts.items() if int(key[0]) == length]
        add(
            f"regional_gc6_{length}nt_exploratory_family_size",
            "PASS" if counts and set(counts) == {expected} else "FAIL",
            expected if counts and set(counts) == {expected} else str(sorted(set(counts))),
            f"families={len(counts)}; expected_each={expected}",
        )
    add(
        "regional_gc9_14_regression",
        "PASS" if regional_regression_pass else "FAIL",
        regional_regression_n,
        f"maximum_absolute_difference={regional_regression_max:.17g}; tolerance={REGRESSION_TOLERANCE}",
    )
    return rows


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        if value.is_integer():
            return str(int(value))
        return format(value, ".15g")
    return str(value)


def write_table(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def read_eligibility(path: Path) -> list[dict[str, str]]:
    return read_tsv(path)


def validate_stage02_references(
    pair_reference: Path,
    across_reference: Path,
) -> dict[str, str]:
    """Validate immutable Stage 02 regression artifacts and record their identities."""
    references = {
        "stage02_pair_reference": pair_reference,
        "stage02_across_reference": across_reference,
    }
    identities: dict[str, str] = {}
    for label, path in references.items():
        if not path.is_file():
            raise Stage07Error(f"required Stage 02 regression reference is missing: {path}")
        digest = hashlib.sha256()
        try:
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
        except OSError as exc:
            raise Stage07Error(f"cannot read Stage 02 regression reference {path}: {exc}") from exc
        identities[f"{label}_path"] = str(path)
        identities[f"{label}_sha256"] = digest.hexdigest()
    return identities


def run(
    legacy_core: Path,
    stage02_pair_reference: Path,
    stage02_across_reference: Path,
    output_root: Path,
) -> tuple[float, bool]:
    started = time.monotonic()
    stage02_identities = validate_stage02_references(
        stage02_pair_reference, stage02_across_reference
    )
    eligibility = read_eligibility(legacy_core / "results/descriptive/eligibility.tsv")
    primary_pairs = {
        (row["sample"], row["analysis_unit"])
        for row in eligibility
        if is_true(row.get("primary_eligible", ""))
    }
    backgrounds = {
        pair: parse_fasta(
            legacy_core
            / "references/consensus"
            / f"{pair[0]}.{pair[1]}.final.background_masked.fa"
        )
        for pair in primary_pairs
    }
    pair_rows, gc_pair_rows, stats = build_pair_tables(
        eligibility,
        iter_feature_rows(legacy_core, {row["sample"] for row in eligibility}),
        backgrounds,
    )
    sample_rows, summary_rows = aggregate_positional(pair_rows)
    gc_sample_rows, gc_summary_rows = aggregate_gc(gc_pair_rows)
    literature_rows = build_literature_validation(summary_rows, gc_summary_rows)
    regional_pair_rows = build_regional_gc_pair_rows(pair_rows)
    regional_sample_rows, regional_summary_rows, regional_discovery_rows = aggregate_regional_gc(
        regional_pair_rows, literature_rows
    )
    regional_regression_pass, regional_regression_max, regional_regression_n = (
        regional_gc9_14_regression(
            regional_pair_rows,
            regional_sample_rows,
            regional_summary_rows,
            gc_pair_rows,
            gc_sample_rows,
            gc_summary_rows,
        )
    )
    discovery_rows = [
        row for row in summary_rows
        if row["strand"] == "antisense"
        and 3 <= int(row["position_5p"]) <= int(row["length"]) - 2
    ]
    sense_rows = [
        row for row in summary_rows
        if row["strand"] == "sense"
        and 3 <= int(row["position_5p"]) <= int(row["length"]) - 2
    ]
    regression_rows, regression_pass, regression_max, regression_n = stage02_terminal_regression(
        pair_rows,
        summary_rows,
        stage02_pair_reference,
        stage02_across_reference,
    )
    qc_rows = build_qc(
        stats,
        regression_pass,
        regression_max,
        regression_n,
        regional_summary_rows,
        regional_discovery_rows,
        regional_regression_pass,
        regional_regression_max,
        regional_regression_n,
    )
    provenance_rows = [
        {"parameter": "stage", "value": "07_empirical_sequence"},
        {"parameter": "target_lengths", "value": "23,24"},
        {"parameter": "strands", "value": "sense,antisense"},
        {"parameter": "weighting_modes", "value": "unique_sequence,abundance,accumulation"},
        {"parameter": "bootstrap_replicates", "value": BOOTSTRAP_REPLICATES},
        {"parameter": "bootstrap_seed", "value": BOOTSTRAP_SEED},
        {"parameter": "bootstrap_method", "value": "sample-clustered percentile"},
        {"parameter": "ci_level", "value": CI_LEVEL},
        {"parameter": "paired_test", "value": "two-sided exact sign test"},
        {"parameter": "literature_family", "value": "12 tests; Benjamini-Hochberg"},
        {"parameter": "discovery_primary_correction", "value": "Benjamini-Yekutieli"},
        {"parameter": "discovery_sensitivity_correction", "value": "Benjamini-Hochberg"},
        {"parameter": "terminal_regression_tolerance", "value": REGRESSION_TOLERANCE},
        {"parameter": "regional_gc_window_width_nt", "value": REGIONAL_GC_WIDTH},
        {"parameter": "regional_gc_primary_correction", "value": "Benjamini-Yekutieli"},
        {"parameter": "regional_gc_sensitivity_correction", "value": "Benjamini-Hochberg"},
        {"parameter": "regional_gc9_14_regression_maximum_absolute_difference", "value": regional_regression_max},
        {"parameter": "legacy_core", "value": str(legacy_core)},
        {"parameter": "stage02_pair_reference_path", "value": stage02_identities["stage02_pair_reference_path"]},
        {"parameter": "stage02_pair_reference_sha256", "value": stage02_identities["stage02_pair_reference_sha256"]},
        {"parameter": "stage02_across_reference_path", "value": stage02_identities["stage02_across_reference_path"]},
        {"parameter": "stage02_across_reference_sha256", "value": stage02_identities["stage02_across_reference_sha256"]},
        {"parameter": "analysis_scope", "value": "representation/accumulation association; not efficacy"},
    ]

    outputs = {
        "positional_by_pair.tsv": (pair_rows, PAIR_FIELDS),
        "positional_by_sample.tsv": (sample_rows, SAMPLE_FIELDS),
        "positional_summary.tsv": (summary_rows, SUMMARY_FIELDS),
        "gc9_14_by_pair.tsv": (gc_pair_rows, GC_PAIR_FIELDS),
        "gc9_14_by_sample.tsv": (gc_sample_rows, GC_SAMPLE_FIELDS),
        "gc9_14_summary.tsv": (gc_summary_rows, GC_SUMMARY_FIELDS),
        "regional_gc6_by_pair.tsv": (regional_pair_rows, REGIONAL_PAIR_FIELDS),
        "regional_gc6_by_sample.tsv": (regional_sample_rows, REGIONAL_SAMPLE_FIELDS),
        "regional_gc6_summary.tsv": (regional_summary_rows, REGIONAL_SUMMARY_FIELDS),
        "regional_gc6_discovery.tsv": (regional_discovery_rows, REGIONAL_SUMMARY_FIELDS),
        "literature_validation.tsv": (literature_rows, LITERATURE_FIELDS),
        "discovery_summary.tsv": (discovery_rows, SUMMARY_FIELDS),
        "sense_comparator.tsv": (sense_rows, SUMMARY_FIELDS),
        "qc/stage07_accounting.tsv": (qc_rows, QC_FIELDS),
        "qc/stage02_terminal_regression.tsv": (regression_rows, REGRESSION_FIELDS),
        "provenance/stage07_manifest.tsv": (provenance_rows, ["parameter", "value"]),
    }
    for relative_path, (rows, fields) in outputs.items():
        write_table(output_root / relative_path, rows, fields)
    failed = any(row["status"] == "FAIL" for row in qc_rows)
    return time.monotonic() - started, failed


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", required=True, type=Path)
    parser.add_argument("--stage02-pair-reference", required=True, type=Path)
    parser.add_argument("--stage02-across-reference", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        elapsed, failed = run(
            args.legacy_core.resolve(),
            args.stage02_pair_reference.resolve(),
            args.stage02_across_reference.resolve(),
            args.output_root.resolve(),
        )
    except (OSError, ValueError, Stage07Error) as exc:
        print(f"Stage 07 failed: {exc}", file=sys.stderr)
        return 1
    print(f"Stage 07 completed in {elapsed:.3f} seconds", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
