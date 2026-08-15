#!/usr/bin/env python3
"""Canonical Stage 02 matched terminal-nucleotide enrichment tables."""

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
from typing import Iterable, Iterator


NUCLEOTIDES = ("A", "C", "G", "T")
STRANDS = ("sense", "antisense")
SCOPES = ("sense", "antisense", "combined")
WEIGHTING_MODES = ("abundance", "unique_sequence")
EXPECTED_CATEGORIES = {
    "mapping_mode": {"exact", "1mm"},
    "virus_assignment": {"assigned", "ambiguous_multi_virus"},
    "strand": {"sense", "antisense", "ambiguous"},
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class Stage02Config:
    target_lengths: tuple[int, ...]
    terminal_positions: tuple[str, ...]
    bootstrap_replicates: int
    random_seed: int
    ci_method: str
    ci_level: float
    frequency_sum_tolerance: float


def load_config(path: Path) -> Stage02Config:
    data = json.loads(path.read_text())["stage02"]
    config = Stage02Config(
        tuple(int(value) for value in data["target_lengths"]),
        tuple(data["terminal_positions"]),
        int(data["bootstrap_replicates"]),
        int(data["random_seed"]),
        str(data["ci_method"]),
        float(data["ci_level"]),
        float(data["frequency_sum_tolerance"]),
    )
    if config.target_lengths != (23, 24):
        raise ValueError("canonical Stage 02 target_lengths must be [23, 24]")
    if config.terminal_positions != ("5p1", "5p2", "3p2", "3p1"):
        raise ValueError("canonical Stage 02 terminal positions differ from the specification")
    if config.bootstrap_replicates <= 0 or config.ci_method != "percentile":
        raise ValueError("invalid canonical Stage 02 bootstrap configuration")
    if not 0 < config.ci_level < 1 or config.frequency_sum_tolerance <= 0:
        raise ValueError("invalid Stage 02 CI level or frequency tolerance")
    return config


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def terminal_bases(sequence: str) -> dict[str, str]:
    if len(sequence) < 2:
        raise ValueError("terminal extraction requires a sequence of length at least two")
    return {"5p1": sequence[0], "5p2": sequence[1], "3p2": sequence[-2], "3p1": sequence[-1]}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def safe_fraction(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator != 0 else None


def enrichment_ratio(observed: float | None, expected: float | None, observed_total: float) -> float | None:
    if observed_total == 0 or observed is None or expected is None or expected == 0:
        return None
    return observed / expected


def finite(values: Iterable[float | None]) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(float(value))]


def median_or_none(values: Iterable[float | None]) -> float | None:
    usable = finite(values)
    return statistics.median(usable) if usable else None


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def sample_clustered_bootstrap(
    values_by_sample: dict[str, list[float]], replicates: int, seed: int, level: float
) -> tuple[float | None, float | None, int]:
    samples = sorted(values_by_sample)
    if not samples:
        return None, None, 0
    rng = random.Random(seed)
    output: list[float] = []
    for _ in range(replicates):
        selected = [rng.choice(samples) for _ in samples]
        sample_medians = [statistics.median(values_by_sample[sample]) for sample in selected if values_by_sample[sample]]
        if sample_medians:
            output.append(float(statistics.median(sample_medians)))
    if not output:
        return None, None, 0
    alpha = (1 - level) / 2
    return percentile(output, alpha), percentile(output, 1 - alpha), len(output)


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = ((start + 1) + end) / 2
        for position in range(start, end):
            ranks[order[position]] = rank
        start = end
    return ranks


def spearman_rho(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    x, y = average_ranks(left), average_ranks(right)
    mean_x, mean_y = statistics.mean(x), statistics.mean(y)
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - mean_x) ** 2 for a in x) * sum((b - mean_y) ** 2 for b in y))
    return numerator / denominator if denominator else None


def parse_fasta(path: Path) -> list[str]:
    records: list[str] = []
    sequence: list[str] = []
    seen_header = False
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
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
                    raise ValueError(f"sequence before FASTA header in {path} line {line_number}")
                sequence.append(line)
    if seen_header:
        records.append("".join(sequence).upper())
    if not records or any(not sequence for sequence in records):
        raise ValueError(f"empty FASTA record in {path}")
    return records


def enumerate_background(records: list[str], length: int) -> dict[str, object]:
    counts = {strand: Counter() for strand in STRANDS}
    valid = 0
    candidate = 0
    for record in records:
        starts = max(0, len(record) - length + 1)
        candidate += starts
        for start in range(starts):
            window = record[start:start + length]
            if any(base not in NUCLEOTIDES for base in window):
                continue
            valid += 1
            for position, base in terminal_bases(window).items():
                counts["sense"][(position, base)] += 1
            antisense = reverse_complement(window)
            for position, base in terminal_bases(antisense).items():
                counts["antisense"][(position, base)] += 1
    return {"counts": counts, "valid": valid, "candidate": candidate, "excluded": candidate - valid}


def iter_feature_rows(root: Path, samples: Iterable[str]) -> Iterator[dict[str, str]]:
    for sample in sorted(samples):
        path = root / "tables" / sample / f"{sample}.read_level_features.tsv.gz"
        with gzip.open(path, "rt", newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def aggregate_enrichment(
    pair_rows: list[dict[str, object]], config: Stage02Config
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        key = tuple(row[field] for field in ("length", "strand_scope", "weighting_mode", "terminal_position", "nucleotide"))
        grouped[key].append(row)
    sample_rows: list[dict[str, object]] = []
    across_rows: list[dict[str, object]] = []
    for key in sorted(grouped):
        length, scope, mode, position, nucleotide = key
        rows = grouped[key]
        by_sample: dict[str, list[float]] = defaultdict(list)
        for row in rows:
            value = row["enrichment_ratio"]
            if value is not None and math.isfinite(float(value)):
                by_sample[str(row["sample"])].append(float(value))
        for sample in sorted({str(row["sample"]) for row in rows}):
            values = by_sample.get(sample, [])
            sample_rows.append({
                "sample": sample, "length": length, "strand_scope": scope,
                "weighting_mode": mode, "terminal_position": position,
                "nucleotide": nucleotide,
                "sample_enrichment_median": statistics.median(values) if values else None,
                "n_virus_units": len(values),
            })
        ci_low, ci_high, valid = sample_clustered_bootstrap(
            dict(by_sample), config.bootstrap_replicates, config.random_seed, config.ci_level
        )
        pair_values = [float(row["enrichment_ratio"]) for row in rows if row["enrichment_ratio"] is not None and math.isfinite(float(row["enrichment_ratio"]))]
        across_rows.append({
            "length": length, "strand_scope": scope, "weighting_mode": mode,
            "terminal_position": position, "nucleotide": nucleotide,
            "sample_balanced_median_enrichment_ratio": median_or_none(statistics.median(values) for values in by_sample.values()),
            "ci_low": ci_low, "ci_high": ci_high, "n_samples": len(by_sample),
            "n_sample_virus_units": len(pair_values),
            "pair_median_enrichment_ratio": median_or_none(pair_values),
            "bootstrap_replicates_requested": config.bootstrap_replicates,
            "bootstrap_replicates_valid": valid, "bootstrap_seed": config.random_seed,
            "ci_method": config.ci_method, "ci_level": config.ci_level,
        })
    return sample_rows, across_rows


def pooled_abundance(pair_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        if row["weighting_mode"] == "abundance":
            key = tuple(row[field] for field in ("length", "strand_scope", "terminal_position", "nucleotide"))
            grouped[key].append(row)
    output: list[dict[str, object]] = []
    for key in sorted(grouped):
        rows = [row for row in grouped[key] if row["expected_fraction"] is not None and float(row["observed_total_weight"]) > 0]
        total = sum(float(row["observed_total_weight"]) for row in rows)
        observed_weight = sum(float(row["observed_terminal_weight"]) for row in rows)
        expected_weight = sum(float(row["observed_total_weight"]) * float(row["expected_fraction"]) for row in rows)
        observed_fraction = safe_fraction(observed_weight, total)
        expected_fraction = safe_fraction(expected_weight, total)
        ratio = None if observed_fraction is None or expected_fraction in (None, 0) else observed_fraction / expected_fraction
        output.append({
            "length": key[0], "strand_scope": key[1], "terminal_position": key[2],
            "nucleotide": key[3], "pooled_abundance_observed_fraction": observed_fraction,
            "pooled_abundance_expected_fraction": expected_fraction,
            "pooled_abundance_enrichment_ratio": ratio,
            "pooled_observed_total_weight": total,
            "n_samples": len({row["sample"] for row in rows}),
            "n_sample_virus_units": len(rows), "analysis_role": "secondary_descriptive",
        })
    return output


def compare_lengths(across_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    lookup = {
        (row["length"], row["strand_scope"], row["weighting_mode"], row["terminal_position"], row["nucleotide"]): row["sample_balanced_median_enrichment_ratio"]
        for row in across_rows
    }
    output = []
    for scope in ("antisense", "combined"):
        for mode in WEIGHTING_MODES:
            left, right = [], []
            for position in ("5p1", "5p2", "3p2", "3p1"):
                for nucleotide in NUCLEOTIDES:
                    value23 = lookup.get((23, scope, mode, position, nucleotide))
                    value24 = lookup.get((24, scope, mode, position, nucleotide))
                    if value23 is not None and value24 is not None:
                        left.append(float(value23)); right.append(float(value24))
            output.append({
                "strand_scope": scope, "weighting_mode": mode,
                "n_matched_features": len(left), "spearman_rho_23_24": spearman_rho(left, right),
                "input_metric": "sample_balanced_median_enrichment_ratio",
            })
    return output


def analyse_stage02(
    eligibility: list[dict[str, str]], feature_rows: Iterable[dict[str, str]],
    backgrounds: dict[tuple[str, str], list[str]], config: Stage02Config,
) -> dict[str, list[dict[str, object]]]:
    samples = sorted({row["sample"] for row in eligibility})
    metadata = {(row["sample"], row["analysis_unit"]): row for row in eligibility if is_true(row["primary_eligible"])}
    pairs = set(metadata)
    abundance_counts = {pair: Counter() for pair in pairs}
    abundance_totals = {pair: Counter() for pair in pairs}
    unique_counts = {pair: Counter() for pair in pairs}
    unique_totals = {pair: Counter() for pair in pairs}
    sample_unique: dict[tuple[str, str, int, str], set[str]] = defaultdict(set)
    current_sample: str | None = None
    completed_samples: set[str] = set()
    rows_examined = retained_rows = length_mismatches = unexpected_bases = 0
    category_values: dict[str, set[str]] = defaultdict(set)

    def add_sequence(counter: Counter, totals: Counter, pair: tuple[str, str], length: int, strand: str, sequence: str, weight: float) -> None:
        totals[pair][(length, strand)] += weight
        for position, base in terminal_bases(sequence).items():
            counter[pair][(length, strand, position, base)] += weight

    def flush_unique() -> None:
        for (sample, unit, length, strand), sequences in sample_unique.items():
            pair = (sample, unit)
            for sequence in sequences:
                add_sequence(unique_counts, unique_totals, pair, length, strand, sequence, 1.0)
        sample_unique.clear()

    for row in feature_rows:
        rows_examined += 1
        row_sample = row.get("sample", "")
        if current_sample is None:
            current_sample = row_sample
        elif row_sample != current_sample:
            flush_unique(); completed_samples.add(current_sample)
            if row_sample in completed_samples:
                raise ValueError("feature rows must be grouped by sample")
            current_sample = row_sample
        for column in EXPECTED_CATEGORIES:
            category_values[column].add(row.get(column, ""))
        pair = (row_sample, row.get("virus", ""))
        if pair not in pairs or row.get("mapping_mode") != "exact" or row.get("virus_assignment") != "assigned" or row.get("strand") not in STRANDS:
            continue
        try:
            length = int(row["length"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid declared length at read-level row {rows_examined}") from exc
        if length not in config.target_lengths:
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
            raise ValueError(f"invalid count at read-level row {rows_examined}") from exc
        if count < 0 or not math.isfinite(count):
            raise ValueError(f"invalid count at read-level row {rows_examined}")
        retained_rows += 1
        strand = row["strand"]
        add_sequence(abundance_counts, abundance_totals, pair, length, strand, sequence, count)
        sample_unique[(pair[0], pair[1], length, strand)].add(sequence)
    flush_unique()

    observed_rows: list[dict[str, object]] = []
    observed_lookup: dict[tuple[object, ...], dict[str, object]] = {}
    for pair in sorted(pairs):
        meta = metadata[pair]
        for length in config.target_lengths:
            for mode, counters, totals in (
                ("abundance", abundance_counts, abundance_totals),
                ("unique_sequence", unique_counts, unique_totals),
            ):
                sense_total = float(totals[pair][(length, "sense")])
                antisense_total = float(totals[pair][(length, "antisense")])
                combined_total = sense_total + antisense_total
                w_sense = safe_fraction(sense_total, combined_total)
                w_antisense = safe_fraction(antisense_total, combined_total)
                for scope in SCOPES:
                    total = combined_total if scope == "combined" else float(totals[pair][(length, scope)])
                    for position in config.terminal_positions:
                        for nucleotide in NUCLEOTIDES:
                            if scope == "combined":
                                weight = float(counters[pair][(length, "sense", position, nucleotide)] + counters[pair][(length, "antisense", position, nucleotide)])
                            else:
                                weight = float(counters[pair][(length, scope, position, nucleotide)])
                            row = {
                                "sample": pair[0], "analysis_unit": pair[1],
                                "biological_virus": meta["biological_virus"], "polarity": meta["polarity"],
                                "length": length, "strand_scope": scope, "weighting_mode": mode,
                                "terminal_position": position, "nucleotide": nucleotide,
                                "observed_terminal_weight": weight, "observed_total_weight": total,
                                "observed_fraction": safe_fraction(weight, total),
                                "observed_strand_weight_sense": w_sense,
                                "observed_strand_weight_antisense": w_antisense,
                            }
                            observed_rows.append(row)
                            observed_lookup[(pair, length, scope, mode, position, nucleotide)] = row

    background_stats = {(pair, length): enumerate_background(backgrounds[pair], length) for pair in pairs for length in config.target_lengths}
    expected_rows: list[dict[str, object]] = []
    expected_lookup: dict[tuple[object, ...], dict[str, object]] = {}
    for pair in sorted(pairs):
        meta = metadata[pair]
        for length in config.target_lengths:
            stats = background_stats[(pair, length)]
            valid = int(stats["valid"])
            counts = stats["counts"]
            for mode in WEIGHTING_MODES:
                combined_observed = observed_lookup[(pair, length, "combined", mode, "5p1", "A")]
                w_sense = combined_observed["observed_strand_weight_sense"]
                w_antisense = combined_observed["observed_strand_weight_antisense"]
                for scope in SCOPES:
                    for position in config.terminal_positions:
                        for nucleotide in NUCLEOTIDES:
                            if valid == 0:
                                fraction = None
                            elif scope == "combined":
                                sense_fraction = counts["sense"][(position, nucleotide)] / valid
                                antisense_fraction = counts["antisense"][(position, nucleotide)] / valid
                                fraction = None if w_sense is None else float(w_sense) * sense_fraction + float(w_antisense) * antisense_fraction
                            else:
                                fraction = counts[scope][(position, nucleotide)] / valid
                            terminal_weight = None if fraction is None else fraction * valid
                            row = {
                                "sample": pair[0], "analysis_unit": pair[1],
                                "biological_virus": meta["biological_virus"], "polarity": meta["polarity"],
                                "length": length, "strand_scope": scope, "weighting_mode": mode,
                                "terminal_position": position, "nucleotide": nucleotide,
                                "background_fasta_records": len(backgrounds[pair]),
                                "candidate_background_windows": stats["candidate"],
                                "valid_background_windows": valid,
                                "excluded_background_windows": stats["excluded"],
                                "expected_terminal_weight": terminal_weight,
                                "expected_total_weight": valid if fraction is not None else None,
                                "expected_fraction": fraction,
                                "observed_strand_weight_sense": w_sense if scope == "combined" else None,
                                "observed_strand_weight_antisense": w_antisense if scope == "combined" else None,
                            }
                            expected_rows.append(row)
                            expected_lookup[(pair, length, scope, mode, position, nucleotide)] = row

    pair_enrichment: list[dict[str, object]] = []
    for observed in observed_rows:
        pair = (str(observed["sample"]), str(observed["analysis_unit"]))
        key = (pair, observed["length"], observed["strand_scope"], observed["weighting_mode"], observed["terminal_position"], observed["nucleotide"])
        expected = expected_lookup[key]
        pair_enrichment.append(observed | {
            "expected_fraction": expected["expected_fraction"],
            "valid_background_windows": expected["valid_background_windows"],
            "enrichment_ratio": enrichment_ratio(observed["observed_fraction"], expected["expected_fraction"], float(observed["observed_total_weight"])),
        })

    sample_enrichment, across_enrichment = aggregate_enrichment(pair_enrichment, config)
    pooled = pooled_abundance(pair_enrichment)
    comparisons = compare_lengths(across_enrichment)

    observed_deviations = []
    expected_deviations = []
    for rows, field, target in ((observed_rows, "observed_fraction", observed_deviations), (expected_rows, "expected_fraction", expected_deviations)):
        grouped: dict[tuple[object, ...], list[float]] = defaultdict(list)
        for row in rows:
            value = row[field]
            if value is not None:
                grouped[tuple(row[key] for key in ("sample", "analysis_unit", "length", "strand_scope", "weighting_mode", "terminal_position"))].append(float(value))
        target.extend(abs(sum(values) - 1) for values in grouped.values())

    qc: list[dict[str, object]] = []
    def q(metric: str, value: object, status: str = "INFO", details: str = "") -> None:
        qc.append({"metric": metric, "status": status, "value": value, "details": details})
    q("samples_represented", len(samples), "PASS")
    q("primary_eligible_sample_virus_units", len(pairs), "PASS")
    q("read_level_rows_examined", rows_examined)
    q("exact_assigned_eligible_23_24_rows_retained", retained_rows)
    for mode, totals in (("abundance", abundance_totals), ("unique_sequence", unique_totals)):
        for length in config.target_lengths:
            for strand in STRANDS:
                q(f"retained_{mode}_{length}nt_{strand}", sum(totals[pair][(length, strand)] for pair in pairs))
    zero_observed = []
    for pair in sorted(pairs):
        for mode, totals in (("abundance", abundance_totals), ("unique_sequence", unique_totals)):
            for length in config.target_lengths:
                for scope in SCOPES:
                    total = sum(totals[pair][(length, strand)] for strand in STRANDS) if scope == "combined" else totals[pair][(length, scope)]
                    if total == 0:
                        zero_observed.append(f"{pair[0]}:{pair[1]}:{length}:{scope}:{mode}")
    q("zero_signal_observed_populations", len(zero_observed), "WARN" if zero_observed else "PASS", ",".join(zero_observed))
    q("declared_length_sequence_length_mismatches", length_mismatches, "FAIL" if length_mismatches else "PASS")
    q("unexpected_observed_bases", unexpected_bases, "FAIL" if unexpected_bases else "PASS")
    for column, expected in EXPECTED_CATEGORIES.items():
        unexpected = sorted(category_values[column] - expected)
        q(f"unexpected_{column}", len(unexpected), "WARN" if unexpected else "PASS", f"observed={sorted(category_values[column])}; unexpected={unexpected}")
    q("background_fasta_records_examined", sum(len(records) for records in backgrounds.values()))
    for length in config.target_lengths:
        q(f"valid_{length}nt_background_windows", sum(int(background_stats[(pair, length)]["valid"]) for pair in pairs))
        q(f"excluded_{length}nt_background_windows", sum(int(background_stats[(pair, length)]["excluded"]) for pair in pairs))
    zero_background = [f"{pair[0]}:{pair[1]}:{length}" for pair in sorted(pairs) for length in config.target_lengths if int(background_stats[(pair, length)]["valid"]) == 0]
    q("units_with_zero_valid_background_windows", len(zero_background), "WARN" if zero_background else "PASS", ",".join(zero_background))
    for pair in sorted(pairs):
        for length in config.target_lengths:
            stats = background_stats[(pair, length)]
            q("background_windows_by_unit", stats["valid"], "INFO", f"sample={pair[0]}; analysis_unit={pair[1]}; length={length}; excluded={stats['excluded']}")
    q("maximum_observed_frequency_sum_deviation", max(observed_deviations, default=0.0), "PASS" if max(observed_deviations, default=0.0) <= config.frequency_sum_tolerance else "FAIL")
    q("maximum_expected_frequency_sum_deviation", max(expected_deviations, default=0.0), "PASS" if max(expected_deviations, default=0.0) <= config.frequency_sum_tolerance else "FAIL")
    finite_enrichment = sum(row["enrichment_ratio"] is not None and math.isfinite(float(row["enrichment_ratio"])) for row in pair_enrichment)
    q("finite_pair_level_enrichment_values", finite_enrichment)
    q("non_finite_pair_level_enrichment_values", len(pair_enrichment) - finite_enrichment)
    sample_counts = [int(row["n_samples"]) for row in across_enrichment]
    q("canonical_feature_contributing_samples_min", min(sample_counts, default=0))
    q("canonical_feature_contributing_samples_max", max(sample_counts, default=0))
    pooled_totals = [float(row["pooled_observed_total_weight"]) for row in pooled]
    q("pooled_abundance_feature_rows", len(pooled))
    q("pooled_observed_total_weight_min", min(pooled_totals, default=0.0))
    q("pooled_observed_total_weight_max", max(pooled_totals, default=0.0))

    return {
        "qc": qc, "observed": observed_rows, "expected": expected_rows,
        "enrichment_pair": pair_enrichment, "enrichment_sample": sample_enrichment,
        "enrichment_across": across_enrichment, "pooled": pooled, "comparisons": comparisons,
    }


FIELDS = {
    "qc": ["metric", "status", "value", "details"],
    "observed": ["sample", "analysis_unit", "biological_virus", "polarity", "length", "strand_scope", "weighting_mode", "terminal_position", "nucleotide", "observed_terminal_weight", "observed_total_weight", "observed_fraction", "observed_strand_weight_sense", "observed_strand_weight_antisense"],
    "expected": ["sample", "analysis_unit", "biological_virus", "polarity", "length", "strand_scope", "weighting_mode", "terminal_position", "nucleotide", "background_fasta_records", "candidate_background_windows", "valid_background_windows", "excluded_background_windows", "expected_terminal_weight", "expected_total_weight", "expected_fraction", "observed_strand_weight_sense", "observed_strand_weight_antisense"],
    "enrichment_pair": ["sample", "analysis_unit", "biological_virus", "polarity", "length", "strand_scope", "weighting_mode", "terminal_position", "nucleotide", "observed_terminal_weight", "observed_total_weight", "observed_fraction", "expected_fraction", "valid_background_windows", "enrichment_ratio"],
    "enrichment_sample": ["sample", "length", "strand_scope", "weighting_mode", "terminal_position", "nucleotide", "sample_enrichment_median", "n_virus_units"],
    "enrichment_across": ["length", "strand_scope", "weighting_mode", "terminal_position", "nucleotide", "sample_balanced_median_enrichment_ratio", "ci_low", "ci_high", "n_samples", "n_sample_virus_units", "pair_median_enrichment_ratio", "bootstrap_replicates_requested", "bootstrap_replicates_valid", "bootstrap_seed", "ci_method", "ci_level"],
    "pooled": ["length", "strand_scope", "terminal_position", "nucleotide", "pooled_abundance_observed_fraction", "pooled_abundance_expected_fraction", "pooled_abundance_enrichment_ratio", "pooled_observed_total_weight", "n_samples", "n_sample_virus_units", "analysis_role"],
    "comparisons": ["strand_scope", "weighting_mode", "n_matched_features", "spearman_rho_23_24", "input_metric"],
}

PATHS = {
    "qc": "qc/stage02_accounting.tsv",
    "observed": "observed/terminal_observed_by_pair.tsv",
    "expected": "background/terminal_expected_by_pair.tsv",
    "enrichment_pair": "enrichment/terminal_enrichment_by_pair.tsv",
    "enrichment_sample": "enrichment/terminal_enrichment_by_sample.tsv",
    "enrichment_across": "enrichment/terminal_enrichment_across_dataset.tsv",
    "pooled": "enrichment/terminal_enrichment_pooled_abundance.tsv",
    "comparisons": "comparisons/enrichment_23_vs_24.tsv",
}


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
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def read_eligibility(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(legacy_core: Path, config_path: Path, output_root: Path) -> tuple[float, bool]:
    started = time.monotonic()
    config = load_config(config_path)
    eligibility = read_eligibility(legacy_core / "results/descriptive/eligibility.tsv")
    primary_pairs = {(row["sample"], row["analysis_unit"]) for row in eligibility if is_true(row["primary_eligible"])}
    backgrounds = {
        pair: parse_fasta(legacy_core / "references/consensus" / f"{pair[0]}.{pair[1]}.final.background_masked.fa")
        for pair in primary_pairs
    }
    results = analyse_stage02(
        eligibility, iter_feature_rows(legacy_core, {row["sample"] for row in eligibility}), backgrounds, config
    )
    for key in PATHS:
        write_table(output_root / PATHS[key], results[key], FIELDS[key])
    failed = any(row["status"] == "FAIL" for row in results["qc"])
    return time.monotonic() - started, failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args(argv)
    elapsed, failed = run(args.legacy_core.resolve(), args.config.resolve(), args.output_root.resolve())
    print(f"Stage 02 completed in {elapsed:.3f} seconds", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
