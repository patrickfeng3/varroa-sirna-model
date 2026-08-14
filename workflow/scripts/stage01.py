#!/usr/bin/env python3
"""Canonical Stage 01 viral length landscape and fixed 23/24 tables."""

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


WEIGHTING_MODES = ("abundance", "unique_sequence")
STRANDS = ("sense", "antisense")
FIXED_METRICS = (
    "sense_fraction_23",
    "antisense_fraction_23",
    "sense_fraction_24",
    "antisense_fraction_24",
    "delta_antisense_fraction_24_minus_23",
    "length23_fraction_among_23_24",
    "length24_fraction_among_23_24",
)
EXPECTED_CATEGORIES = {
    "mapping_mode": {"exact", "1mm"},
    "virus_assignment": {"assigned", "ambiguous_multi_virus"},
    "strand": {"sense", "antisense", "ambiguous"},
}


@dataclass(frozen=True)
class Stage01Config:
    length_min: int
    length_max: int
    bootstrap_replicates: int
    random_seed: int
    ci_method: str
    ci_level: float

    @property
    def lengths(self) -> range:
        return range(self.length_min, self.length_max + 1)


def load_config(path: Path) -> Stage01Config:
    data = json.loads(path.read_text())["stage01"]
    config = Stage01Config(
        length_min=int(data["length_min"]),
        length_max=int(data["length_max"]),
        bootstrap_replicates=int(data["bootstrap_replicates"]),
        random_seed=int(data["random_seed"]),
        ci_method=str(data["ci_method"]),
        ci_level=float(data["ci_level"]),
    )
    if config.length_min > config.length_max:
        raise ValueError("stage01 length_min exceeds length_max")
    if config.bootstrap_replicates <= 0:
        raise ValueError("bootstrap_replicates must be positive")
    if config.ci_method != "percentile" or not 0 < config.ci_level < 1:
        raise ValueError("Stage 01 requires percentile CI with level between zero and one")
    return config


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def finite(values: Iterable[float | None]) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(float(value))]


def median_or_none(values: Iterable[float | None]) -> float | None:
    usable = finite(values)
    return statistics.median(usable) if usable else None


def safe_fraction(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator != 0 else None


def competition_ranks(counts: dict[int, float]) -> dict[int, int]:
    """Descending standard competition ranks: 1, 2, 2, 4."""
    ordered = sorted(counts.values(), reverse=True)
    first_rank: dict[float, int] = {}
    for index, value in enumerate(ordered, 1):
        first_rank.setdefault(value, index)
    return {length: first_rank[value] for length, value in counts.items()}


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("percentile requires at least one value")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def clustered_bootstrap_ci(
    sample_values: dict[str, float | None], replicates: int, seed: int, level: float
) -> tuple[float | None, float | None, int]:
    usable = {sample: value for sample, value in sample_values.items() if value is not None}
    samples = sorted(sample_values)
    if not usable or not samples:
        return None, None, 0
    rng = random.Random(seed)
    statistics_out: list[float] = []
    for _ in range(replicates):
        selected = [rng.choice(samples) for _ in samples]
        replicate_values = [float(usable[sample]) for sample in selected if sample in usable]
        if replicate_values:
            statistics_out.append(float(statistics.median(replicate_values)))
    alpha = (1 - level) / 2
    if not statistics_out:
        return None, None, 0
    return (
        percentile(statistics_out, alpha),
        percentile(statistics_out, 1 - alpha),
        len(statistics_out),
    )


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


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def iter_feature_rows(root: Path, samples: Iterable[str]) -> Iterator[dict[str, str]]:
    for sample in sorted(samples):
        path = root / "tables" / sample / f"{sample}.read_level_features.tsv.gz"
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                yield row


def analyse_stage01(
    eligibility: list[dict[str, str]],
    feature_rows: Iterable[dict[str, str]],
    config: Stage01Config,
) -> dict[str, list[dict[str, object]]]:
    samples = sorted({row["sample"] for row in eligibility})
    primary_metadata = {
        (row["sample"], row["analysis_unit"]): row
        for row in eligibility if is_true(row["primary_eligible"])
    }
    primary_pairs = set(primary_metadata)
    abundance: dict[tuple[str, str], Counter[tuple[int, str]]] = {
        pair: Counter() for pair in primary_pairs
    }
    unique: dict[tuple[str, str], Counter[tuple[int, str]]] = {
        pair: Counter() for pair in primary_pairs
    }
    sample_unique: dict[tuple[str, str, int, str], set[str]] = defaultdict(set)
    current_sample: str | None = None
    completed_samples: set[str] = set()
    distinct_unique = 0
    unique23 = 0
    unique24 = 0

    def flush_sample_unique() -> None:
        nonlocal distinct_unique, unique23, unique24
        for (sample, unit, length, strand), sequences in sample_unique.items():
            count = len(sequences)
            unique[(sample, unit)][(length, strand)] += count
            distinct_unique += count
            if length == 23:
                unique23 += count
            elif length == 24:
                unique24 += count
        sample_unique.clear()

    category_values: dict[str, set[str]] = defaultdict(set)
    rows_examined = 0
    retained_rows = 0
    outside_rows = 0
    retained_abundance = 0.0

    for row in feature_rows:
        rows_examined += 1
        row_sample = row.get("sample", "")
        if current_sample is None:
            current_sample = row_sample
        elif row_sample != current_sample:
            flush_sample_unique()
            completed_samples.add(current_sample)
            if row_sample in completed_samples:
                raise ValueError("feature rows must be grouped by sample for bounded-memory deduplication")
            current_sample = row_sample
        for column in EXPECTED_CATEGORIES:
            category_values[column].add(row.get(column, ""))
        pair = (row.get("sample", ""), row.get("virus", ""))
        if pair not in primary_pairs:
            continue
        if row.get("mapping_mode") != "exact" or row.get("virus_assignment") != "assigned":
            continue
        if row.get("strand") not in STRANDS:
            continue
        retained_rows += 1
        try:
            length = int(row["length"])
            count = float(row["count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid numeric Stage 01 row {rows_examined}: {exc}") from exc
        if count < 0 or not math.isfinite(count):
            raise ValueError(f"invalid count in Stage 01 row {rows_examined}: {count}")
        if length < config.length_min or length > config.length_max:
            outside_rows += 1
            continue
        retained_abundance += count
        strand = row["strand"]
        abundance[pair][(length, strand)] += count
        sample_unique[(pair[0], pair[1], length, strand)].add(row["sequence"])
    flush_sample_unique()

    weighted = {"abundance": abundance, "unique_sequence": unique}
    length_pair: list[dict[str, object]] = []
    zero_units: dict[str, list[str]] = defaultdict(list)
    for pair in sorted(primary_pairs):
        metadata = primary_metadata[pair]
        for mode in WEIGHTING_MODES:
            counts = {
                length: float(sum(weighted[mode][pair][(length, strand)] for strand in STRANDS))
                for length in config.lengths
            }
            denominator = sum(counts.values())
            ranks = competition_ranks(counts) if denominator else {}
            if denominator == 0:
                zero_units[mode].append(f"{pair[0]}:{pair[1]}")
            for length in config.lengths:
                rank = ranks.get(length)
                length_pair.append({
                    "sample": pair[0], "analysis_unit": pair[1],
                    "biological_virus": metadata["biological_virus"],
                    "polarity": metadata["polarity"], "weighting_mode": mode,
                    "length": length, "length_count": counts[length],
                    "length_fraction": safe_fraction(counts[length], denominator),
                    "length_rank": rank,
                    "top1_indicator": int(rank <= 1) if rank is not None else None,
                    "top3_indicator": int(rank <= 3) if rank is not None else None,
                })

    length_sample: list[dict[str, object]] = []
    for sample in samples:
        for mode in WEIGHTING_MODES:
            for length in config.lengths:
                rows = [
                    row for row in length_pair
                    if row["sample"] == sample and row["weighting_mode"] == mode
                    and row["length"] == length and row["length_fraction"] is not None
                ]
                length_sample.append({
                    "sample": sample, "weighting_mode": mode, "length": length,
                    "median_length_fraction": median_or_none(row["length_fraction"] for row in rows),
                    "median_length_rank": median_or_none(row["length_rank"] for row in rows),
                    "n_virus_units": len(rows),
                })

    length_across: list[dict[str, object]] = []
    for mode in WEIGHTING_MODES:
        for length in config.lengths:
            sample_rows = [
                row for row in length_sample
                if row["weighting_mode"] == mode and row["length"] == length
                and row["median_length_fraction"] is not None
            ]
            sample_values = {row["sample"]: row["median_length_fraction"] for row in sample_rows}
            ci_low, ci_high, valid = clustered_bootstrap_ci(
                sample_values, config.bootstrap_replicates, config.random_seed, config.ci_level
            )
            pair_rows = [
                row for row in length_pair
                if row["weighting_mode"] == mode and row["length"] == length
                and row["length_fraction"] is not None
            ]
            length_across.append({
                "weighting_mode": mode, "length": length,
                "sample_balanced_median_length_fraction": median_or_none(sample_values.values()),
                "ci95_low": ci_low, "ci95_high": ci_high,
                "median_sample_level_rank": median_or_none(row["median_length_rank"] for row in sample_rows),
                "n_samples": len(sample_rows), "n_sample_virus_units": len(pair_rows),
                "pair_top1_frequency": safe_fraction(sum(int(row["top1_indicator"]) for row in pair_rows), len(pair_rows)),
                "pair_top3_frequency": safe_fraction(sum(int(row["top3_indicator"]) for row in pair_rows), len(pair_rows)),
                "bootstrap_replicates_requested": config.bootstrap_replicates,
                "bootstrap_replicates_valid": valid, "random_seed": config.random_seed,
                "ci_method": config.ci_method, "ci_level": config.ci_level,
            })

    counts_pair: list[dict[str, object]] = []
    fractions_pair: list[dict[str, object]] = []
    for pair in sorted(primary_pairs):
        metadata = primary_metadata[pair]
        for mode in WEIGHTING_MODES:
            counter = weighted[mode][pair]
            n23s = float(counter[(23, "sense")]); n23a = float(counter[(23, "antisense")])
            n24s = float(counter[(24, "sense")]); n24a = float(counter[(24, "antisense")])
            n23 = n23s + n23a; n24 = n24s + n24a
            base = {
                "sample": pair[0], "analysis_unit": pair[1],
                "biological_virus": metadata["biological_virus"],
                "polarity": metadata["polarity"], "weighting_mode": mode,
            }
            counts_pair.append(base | {
                "n23_sense": n23s, "n23_antisense": n23a, "n23_total": n23,
                "n24_sense": n24s, "n24_antisense": n24a, "n24_total": n24,
            })
            sense23 = safe_fraction(n23s, n23); antisense23 = safe_fraction(n23a, n23)
            sense24 = safe_fraction(n24s, n24); antisense24 = safe_fraction(n24a, n24)
            delta = antisense24 - antisense23 if antisense23 is not None and antisense24 is not None else None
            total = n23 + n24
            fractions_pair.append(base | {
                "sense_fraction_23": sense23, "antisense_fraction_23": antisense23,
                "sense_fraction_24": sense24, "antisense_fraction_24": antisense24,
                "delta_antisense_fraction_24_minus_23": delta,
                "length23_fraction_among_23_24": safe_fraction(n23, total),
                "length24_fraction_among_23_24": safe_fraction(n24, total),
            })

    fixed_sample: list[dict[str, object]] = []
    for sample in samples:
        for mode in WEIGHTING_MODES:
            rows = [row for row in fractions_pair if row["sample"] == sample and row["weighting_mode"] == mode]
            for metric in FIXED_METRICS:
                values = finite(row[metric] for row in rows)
                fixed_sample.append({
                    "sample": sample, "weighting_mode": mode, "metric": metric,
                    "median_value": statistics.median(values) if values else None,
                    "n_virus_units": len(values),
                })

    fixed_across: list[dict[str, object]] = []
    for mode in WEIGHTING_MODES:
        for metric in FIXED_METRICS:
            sample_rows = [
                row for row in fixed_sample
                if row["weighting_mode"] == mode and row["metric"] == metric
                and row["median_value"] is not None
            ]
            sample_values = {row["sample"]: row["median_value"] for row in sample_rows}
            pair_values = [
                row[metric] for row in fractions_pair
                if row["weighting_mode"] == mode and row[metric] is not None
            ]
            ci_low, ci_high, valid = clustered_bootstrap_ci(
                sample_values, config.bootstrap_replicates, config.random_seed, config.ci_level
            )
            fixed_across.append({
                "weighting_mode": mode, "metric": metric,
                "sample_balanced_median": median_or_none(sample_values.values()),
                "ci95_low": ci_low, "ci95_high": ci_high,
                "pair_balanced_median": median_or_none(pair_values),
                "n_samples": len(sample_values), "n_sample_virus_units": len(pair_values),
                "bootstrap_replicates_requested": config.bootstrap_replicates,
                "bootstrap_replicates_valid": valid, "random_seed": config.random_seed,
                "ci_method": config.ci_method, "ci_level": config.ci_level,
            })

    abundance23 = sum(abundance[pair][(23, strand)] for pair in primary_pairs for strand in STRANDS)
    abundance24 = sum(abundance[pair][(24, strand)] for pair in primary_pairs for strand in STRANDS)
    qc: list[dict[str, object]] = []
    def qc_row(metric: str, value: object, status: str = "INFO", details: str = "") -> None:
        qc.append({"metric": metric, "status": status, "value": value, "details": details})
    qc_row("samples_represented", len(samples), "PASS")
    qc_row("primary_eligible_sample_virus_units", len(primary_pairs), "PASS")
    qc_row("read_level_rows_examined", rows_examined)
    qc_row("exact_assigned_eligible_rows_retained", retained_rows)
    qc_row("abundance_retained_15_35_nt", retained_abundance)
    qc_row("distinct_stage01_sequences_retained", distinct_unique)
    qc_row("rows_outside_15_35_nt", outside_rows)
    qc_row("total_23nt_abundance", abundance23)
    qc_row("total_23nt_unique_sequences", unique23)
    qc_row("total_24nt_abundance", abundance24)
    qc_row("total_24nt_unique_sequences", unique24)
    for mode in WEIGHTING_MODES:
        units = zero_units[mode]
        qc_row(
            f"zero_signal_units_{mode}", len(units), "WARN" if units else "PASS",
            ",".join(units),
        )
    for column, expected in EXPECTED_CATEGORIES.items():
        observed = category_values[column]
        unexpected = sorted(observed - expected)
        qc_row(
            f"unexpected_{column}", len(unexpected), "WARN" if unexpected else "PASS",
            f"observed={sorted(observed)}; unexpected={unexpected}",
        )

    return {
        "qc": qc, "length_pair": length_pair, "length_sample": length_sample,
        "length_across": length_across, "counts_pair": counts_pair,
        "fractions_pair": fractions_pair, "fixed_sample": fixed_sample,
        "fixed_across": fixed_across,
    }


OUTPUTS = {
    "qc": ("qc/stage01_accounting.tsv", ["metric", "status", "value", "details"]),
    "length_pair": ("length_spectrum/length_distribution_by_pair.tsv", ["sample", "analysis_unit", "biological_virus", "polarity", "weighting_mode", "length", "length_count", "length_fraction", "length_rank", "top1_indicator", "top3_indicator"]),
    "length_sample": ("length_spectrum/length_distribution_by_sample.tsv", ["sample", "weighting_mode", "length", "median_length_fraction", "median_length_rank", "n_virus_units"]),
    "length_across": ("length_spectrum/length_distribution_across_dataset.tsv", ["weighting_mode", "length", "sample_balanced_median_length_fraction", "ci95_low", "ci95_high", "median_sample_level_rank", "n_samples", "n_sample_virus_units", "pair_top1_frequency", "pair_top3_frequency", "bootstrap_replicates_requested", "bootstrap_replicates_valid", "random_seed", "ci_method", "ci_level"]),
    "counts_pair": ("fixed_23_24/23_24_counts_by_pair.tsv", ["sample", "analysis_unit", "biological_virus", "polarity", "weighting_mode", "n23_sense", "n23_antisense", "n23_total", "n24_sense", "n24_antisense", "n24_total"]),
    "fractions_pair": ("fixed_23_24/23_24_fractions_by_pair.tsv", ["sample", "analysis_unit", "biological_virus", "polarity", "weighting_mode", *FIXED_METRICS]),
    "fixed_sample": ("fixed_23_24/23_24_by_sample.tsv", ["sample", "weighting_mode", "metric", "median_value", "n_virus_units"]),
    "fixed_across": ("fixed_23_24/23_24_across_dataset.tsv", ["weighting_mode", "metric", "sample_balanced_median", "ci95_low", "ci95_high", "pair_balanced_median", "n_samples", "n_sample_virus_units", "bootstrap_replicates_requested", "bootstrap_replicates_valid", "random_seed", "ci_method", "ci_level"]),
}


def write_table(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def run(legacy_core: Path, config_path: Path, output_root: Path) -> float:
    started = time.monotonic()
    config = load_config(config_path)
    _, eligibility = read_tsv(legacy_core / "results/descriptive/eligibility.tsv")
    samples = {row["sample"] for row in eligibility}
    results = analyse_stage01(eligibility, iter_feature_rows(legacy_core, samples), config)
    for key, (relative, fields) in OUTPUTS.items():
        write_table(output_root / relative, results[key], fields)
    return time.monotonic() - started


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args(argv)
    elapsed = run(args.legacy_core.resolve(), args.config.resolve(), args.output_root.resolve())
    print(f"Stage 01 completed in {elapsed:.3f} seconds", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
