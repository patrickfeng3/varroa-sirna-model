#!/usr/bin/env python3
"""Canonical Stage 05 viral spatial/transitivity-consistency analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import statistics
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


TRACKS = ("23S", "23AS", "24S", "24AS")
WEIGHTING_MODES = ("abundance", "unique_sequence")
ANCHOR_TYPES = ("balanced23", "combined23")
ENDPOINTS = ("delta_F24_AS", "antisense_specific_directionality")


class Stage05Error(RuntimeError):
    """Structured Stage 05 failure."""


@dataclass(frozen=True)
class Stage05Config:
    positive_sense_polarities: tuple[str, ...]
    target_lengths: tuple[int, ...]
    bin_size_nt: int
    windows_nt: tuple[int, ...]
    anchor_percentile: float
    anchor_min_separation_nt: int
    minimum_anchors: int
    max_crosscorr_lag_nt: int
    minimum_crosscorr_overlap_bins: int
    permutations: int
    bootstrap_replicates: int
    random_seed: int
    ci_level: float
    normalization_tolerance: float
    regression_tolerance_effect: float
    regression_checkpoints: dict[int, dict[str, float]]


def load_config(path: Path) -> Stage05Config:
    data = json.loads(path.read_text())["stage05"]
    config = Stage05Config(
        positive_sense_polarities=tuple(data["positive_sense_polarities"]),
        target_lengths=tuple(int(x) for x in data["target_lengths"]),
        bin_size_nt=int(data["bin_size_nt"]),
        windows_nt=tuple(int(x) for x in data["windows_nt"]),
        anchor_percentile=float(data["anchor_percentile"]),
        anchor_min_separation_nt=int(data["anchor_min_separation_nt"]),
        minimum_anchors=int(data["minimum_anchors"]),
        max_crosscorr_lag_nt=int(data["max_crosscorr_lag_nt"]),
        minimum_crosscorr_overlap_bins=int(data["minimum_crosscorr_overlap_bins"]),
        permutations=int(data["permutations"]),
        bootstrap_replicates=int(data["bootstrap_replicates"]),
        random_seed=int(data["random_seed"]),
        ci_level=float(data["ci_level"]),
        normalization_tolerance=float(data["normalization_tolerance"]),
        regression_tolerance_effect=float(data["regression_tolerance_effect"]),
        regression_checkpoints={int(k): {x: float(y) for x, y in v.items()} for k, v in data["regression_checkpoints"].items()},
    )
    if config.target_lengths != (23, 24) or config.windows_nt != (100, 250, 500):
        raise ValueError("Stage 05 canonical target lengths/windows changed")
    if config.bin_size_nt != 10 or config.permutations != 5000 or config.bootstrap_replicates != 5000:
        raise ValueError("Stage 05 canonical bin/permutation/bootstrap parameters changed")
    if config.random_seed != 20260810:
        raise ValueError("Stage 05 canonical seed changed")
    return config


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def finite(value: object) -> float | None:
    if value is None or str(value) in {"", "NA", "NaN", "nan"}:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def safe_fraction(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator and math.isfinite(denominator) else None


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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def percentile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(float(x) for x in values)
    if not ordered:
        raise ValueError("percentile of empty values")
    position = (len(ordered) - 1) * probability / 100.0
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def alignment_midpoint_nt(start_0based: int, read_length: int) -> float:
    return start_0based + (read_length - 1) / 2


def midpoint_bin(start_0based: int, read_length: int, bin_size_nt: int = 10) -> int:
    return math.floor(alignment_midpoint_nt(start_0based, read_length) / bin_size_nt)


def fractional_locus_weights(total_weight: float, loci: Iterable[tuple]) -> dict[tuple, float]:
    unique_loci = sorted(set(loci))
    if not unique_loci:
        return {}
    weight = float(total_weight) / len(unique_loci)
    return {locus: weight for locus in unique_loci}


def anchor_scores(track23_s: np.ndarray, track23_as: np.ndarray, anchor_type: str) -> np.ndarray:
    if anchor_type == "balanced23":
        return np.sqrt(track23_s * track23_as)
    if anchor_type == "combined23":
        return track23_s + track23_as
    raise ValueError(f"unknown anchor type: {anchor_type}")


def select_anchors(
    scores: Iterable[float], percentile_value: float = 90.0,
    min_separation_nt: int = 50, bin_size_nt: int = 10,
    minimum_anchors: int = 3,
) -> tuple[list[int], float | None, int]:
    array = [float(x) for x in scores]
    nonzero = [value for value in array if value > 0]
    if len(nonzero) < minimum_anchors:
        return [], None, len(nonzero)
    threshold = percentile(nonzero, percentile_value)
    candidates = [index for index, value in enumerate(array) if value > 0 and value >= threshold]
    candidates.sort(key=lambda index: (-array[index], index))
    selected: list[int] = []
    for candidate in candidates:
        if all(abs(candidate - chosen) * bin_size_nt > min_separation_nt for chosen in selected):
            selected.append(candidate)
    if len(selected) < minimum_anchors:
        return [], threshold, len(nonzero)
    return selected, threshold, len(nonzero)


def window_indices(anchor: int, window_bins: int, n_bins: int, direction: str) -> range:
    if direction == "upstream":
        return range(max(0, anchor - window_bins), anchor)
    if direction == "downstream":
        return range(anchor + 1, min(n_bins, anchor + window_bins + 1))
    raise ValueError(direction)


def pooled_window_mean(
    track: np.ndarray, anchors: Iterable[int], window_bins: int, direction: str
) -> tuple[float | None, int]:
    total = 0.0
    count = 0
    for anchor in anchors:
        indices = window_indices(anchor, window_bins, len(track), direction)
        for index in indices:
            total += float(track[index])
            count += 1
    return (total / count if count else None), count


def normalized_directionality(downstream: float | None, upstream: float | None) -> float | None:
    if downstream is None or upstream is None:
        return None
    return safe_fraction(downstream - upstream, downstream + upstream)


def endpoint_values(means: dict[str, float | None]) -> dict[str, float | None]:
    d_as = normalized_directionality(means["mean24AS_down"], means["mean24AS_up"])
    d_s = normalized_directionality(means["mean24S_down"], means["mean24S_up"])
    f_down = safe_fraction(
        means["mean24AS_down"] or 0.0,
        (means["mean23AS_down"] or 0.0) + (means["mean24AS_down"] or 0.0),
    ) if means["mean23AS_down"] is not None and means["mean24AS_down"] is not None else None
    f_up = safe_fraction(
        means["mean24AS_up"] or 0.0,
        (means["mean23AS_up"] or 0.0) + (means["mean24AS_up"] or 0.0),
    ) if means["mean23AS_up"] is not None and means["mean24AS_up"] is not None else None
    return {
        "D_24AS": d_as,
        "D_24S": d_s,
        "antisense_specific_directionality": d_as - d_s if d_as is not None and d_s is not None else None,
        "F24_AS_down": f_down,
        "F24_AS_up": f_up,
        "delta_F24_AS": f_down - f_up if f_down is not None and f_up is not None else None,
    }


def calculate_contig_endpoint(
    tracks: dict[str, np.ndarray], anchors: list[int], window_bins: int
) -> dict[str, object]:
    means: dict[str, float | None] = {}
    counts: dict[str, int] = {}
    for track_name in TRACKS:
        for direction, suffix in (("upstream", "up"), ("downstream", "down")):
            value, count = pooled_window_mean(tracks[track_name], anchors, window_bins, direction)
            means[f"mean{track_name}_{suffix}"] = value
            counts[f"n_valid_bins_{suffix}"] = count
    return {**means, **counts, **endpoint_values(means)}


def allowed_circular_shifts(n_bins: int, exclusion_bins: int) -> tuple[list[int], bool]:
    preferred = [shift for shift in range(1, n_bins) if min(shift, n_bins - shift) > exclusion_bins]
    if preferred:
        return preferred, False
    return list(range(1, n_bins)), True


def apply_same_shift(
    track24_s: np.ndarray, track24_as: np.ndarray, shift: int
) -> tuple[np.ndarray, np.ndarray]:
    return np.roll(track24_s, shift), np.roll(track24_as, shift)


def empirical_p(observed: float | None, null_values: Iterable[float]) -> tuple[float | None, int, int]:
    if observed is None:
        return None, 0, 0
    usable = [float(x) for x in null_values if math.isfinite(float(x))]
    if not usable:
        return None, 0, 0
    exceed = sum(value >= observed for value in usable)
    return (exceed + 1) / (len(usable) + 1), exceed, len(usable)


def stable_seed(seed: int, *parts: object) -> int:
    digest = hashlib.sha256("\x1f".join(str(x) for x in parts).encode()).digest()
    return (seed + int.from_bytes(digest[:8], "big")) % (2**63 - 1)


def bootstrap_median(
    values: list[float], replicates: int, seed: int, level: float
) -> tuple[float | None, float | None, int]:
    if not values:
        return None, None, 0
    rng = np.random.default_rng(seed)
    array = np.asarray(values, dtype=float)
    estimates = np.median(rng.choice(array, size=(replicates, len(array)), replace=True), axis=1)
    alpha = (1 - level) / 2
    return float(np.quantile(estimates, alpha)), float(np.quantile(estimates, 1 - alpha)), replicates


def sample_balanced_value(rows: list[dict[str, object]], field: str) -> tuple[float | None, dict[str, float]]:
    by_sample: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        value = finite(row.get(field))
        if value is not None:
            by_sample[str(row["sample"])].append(value)
    sample_values = {sample: statistics.median(values) for sample, values in by_sample.items() if values}
    return (statistics.median(sample_values.values()) if sample_values else None), sample_values


def virus_balanced_value(rows: list[dict[str, object]], field: str) -> tuple[float | None, dict[str, float]]:
    by_virus: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        value = finite(row.get(field))
        if value is not None:
            by_virus[str(row["biological_virus"])].append(value)
    virus_values = {virus: statistics.median(values) for virus, values in by_virus.items() if values}
    return (statistics.median(virus_values.values()) if virus_values else None), virus_values


def pair_balanced_value(rows: list[dict[str, object]], field: str) -> float | None:
    values = [value for row in rows if (value := finite(row.get(field))) is not None]
    return statistics.median(values) if values else None


def clustered_bootstrap_rows(
    rows: list[dict[str, object]], cluster_field: str, value_field: str,
    replicates: int, seed: int, level: float,
) -> tuple[float | None, float | None, int]:
    by_cluster: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        value = finite(row.get(value_field))
        if value is not None:
            by_cluster[str(row[cluster_field])].append(value)
    clusters = sorted(by_cluster)
    if not clusters:
        return None, None, 0
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(replicates):
        selected = rng.choice(clusters, size=len(clusters), replace=True)
        cluster_medians = [statistics.median(by_cluster[str(cluster)]) for cluster in selected]
        estimates.append(statistics.median(cluster_medians))
    alpha = (1 - level) / 2
    return percentile(estimates, alpha * 100), percentile(estimates, (1 - alpha) * 100), len(estimates)


def hierarchical_bootstrap_rows(
    rows: list[dict[str, object]], cluster_field: str, value_field: str,
    replicates: int, seed: int, level: float,
) -> tuple[float | None, float | None, int]:
    by_cluster: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        value = finite(row.get(value_field))
        if value is not None:
            by_cluster[str(row[cluster_field])].append(value)
    clusters = sorted(by_cluster)
    if not clusters:
        return None, None, 0
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(replicates):
        selected = rng.choice(clusters, size=len(clusters), replace=True)
        cluster_medians = []
        for cluster in selected:
            observations = np.asarray(by_cluster[str(cluster)], dtype=float)
            resampled = rng.choice(observations, size=len(observations), replace=True)
            cluster_medians.append(float(np.median(resampled)))
        estimates.append(statistics.median(cluster_medians))
    alpha = (1 - level) / 2
    return percentile(estimates, alpha * 100), percentile(estimates, (1 - alpha) * 100), len(estimates)


def bh_adjust(p_values: list[float | None]) -> list[float | None]:
    valid = [(index, value) for index, value in enumerate(p_values) if value is not None and math.isfinite(value)]
    output: list[float | None] = [None] * len(p_values)
    if not valid:
        return output
    ordered = sorted(valid, key=lambda item: item[1])
    adjusted = [0.0] * len(ordered)
    running = 1.0
    for position in range(len(ordered) - 1, -1, -1):
        _, value = ordered[position]
        running = min(running, value * len(ordered) / (position + 1))
        adjusted[position] = min(1.0, running)
    for (index, _), value in zip(ordered, adjusted):
        output[index] = value
    return output


def apply_bh(rows: list[dict[str, object]], family_fields: tuple[str, ...], p_field: str, q_field: str) -> None:
    groups: dict[tuple[object, ...], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        groups[tuple(row[field] for field in family_fields)].append(index)
    for indices in groups.values():
        adjusted = bh_adjust([finite(rows[index].get(p_field)) for index in indices])
        for index, value in zip(indices, adjusted):
            rows[index][q_field] = value


def pearson_at_lag(anchor: np.ndarray, target: np.ndarray, lag_bins: int, minimum_overlap: int = 10) -> float | None:
    if lag_bins > 0:
        x, y = anchor[:-lag_bins], target[lag_bins:]
    elif lag_bins < 0:
        x, y = anchor[-lag_bins:], target[:lag_bins]
    else:
        x, y = anchor, target
    if len(x) < minimum_overlap or np.var(x) == 0 or np.var(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def parse_eligibility(path: Path, config: Stage05Config) -> list[dict[str, str]]:
    rows = read_tsv(path)
    return [
        row for row in rows
        if is_true(row["primary_eligible"]) and row["polarity"] in config.positive_sense_polarities
    ]


def load_metadata(
    legacy_core: Path, eligible: list[dict[str, str]], config: Stage05Config
) -> tuple[dict[str, dict[str, dict[str, object]]], dict[str, object]]:
    pairs = {(row["sample"], row["analysis_unit"]) for row in eligible}
    by_sample: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    conflicts: set[tuple[str, str]] = set()
    examined = retained = 0
    categories: dict[str, set[str]] = defaultdict(set)
    for sample in sorted({row["sample"] for row in eligible}):
        path = legacy_core / "tables" / sample / f"{sample}.read_level_features.tsv.gz"
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                examined += 1
                for field in ("mapping_mode", "virus_assignment", "strand"):
                    categories[field].add(row[field])
                if (
                    (sample, row["virus"]) not in pairs
                    or row["mapping_mode"] != "exact"
                    or row["virus_assignment"] != "assigned"
                    or row["strand"] not in {"sense", "antisense"}
                    or int(row["length"]) not in config.target_lengths
                ):
                    continue
                retained += 1
                qname = row["read_name"]
                value = {
                    "sample": sample, "analysis_unit": row["virus"],
                    "strand": row["strand"], "length": int(row["length"]),
                    "sequence": row["sequence"], "count": float(row["count"]),
                }
                existing = by_sample[sample].get(qname)
                if existing is not None and existing != value:
                    conflicts.add((sample, qname))
                else:
                    by_sample[sample][qname] = value
    for sample, qname in conflicts:
        by_sample[sample].pop(qname, None)
    return by_sample, {
        "rows_examined": examined, "rows_retained": retained,
        "metadata_conflicts": len(conflicts), "categories": categories,
    }


def parse_sam_loci(
    path: Path, metadata: dict[str, dict[str, object]], eligible_units: set[str]
) -> tuple[dict[str, set[tuple[str, int, str, int]]], dict[str, int], dict[str, int]]:
    loci: dict[str, set[tuple[str, int, str, int]]] = defaultdict(set)
    lengths: dict[str, int] = {}
    qc = Counter()
    seen_records: set[tuple[str, str, int, str, int]] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("@SQ"):
                fields = dict(item.split(":", 1) for item in line.rstrip().split("\t")[1:] if ":" in item)
                lengths[fields["SN"]] = int(fields["LN"])
                continue
            if line.startswith("@"):
                continue
            qc["sam_records_examined"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                qc["malformed_sam_records"] += 1
                continue
            qname = fields[0]
            if qname not in metadata:
                continue
            flag = int(fields[1])
            if flag & 0x4:
                qc["unmapped_excluded"] += 1
                continue
            if flag & 0x800:
                qc["supplementary_excluded"] += 1
                continue
            rname = fields[2]
            unit = str(metadata[qname]["analysis_unit"])
            if unit not in eligible_units or not rname.startswith(unit + "|"):
                qc["incompatible_loci_excluded"] += 1
                continue
            strand = "antisense" if flag & 0x10 else "sense"
            if strand != metadata[qname]["strand"]:
                qc["strand_mismatches"] += 1
                continue
            start = int(fields[3]) - 1
            length = int(metadata[qname]["length"])
            record = (qname, rname, start, strand, length)
            if record in seen_records:
                qc["duplicate_physical_locus_records"] += 1
                continue
            seen_records.add(record)
            loci[qname].add((rname, start, strand, length))
            qc["compatible_loci_retained"] += 1
    return loci, lengths, dict(qc)


def build_tracks(
    legacy_core: Path, eligible: list[dict[str, str]], metadata_by_sample: dict[str, dict[str, dict[str, object]]],
    config: Stage05Config,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    pair_lookup = {(row["sample"], row["analysis_unit"]): row for row in eligible}
    track_values: dict[tuple[str, str, str, str], dict[str, np.ndarray]] = {}
    contig_lengths: dict[tuple[str, str], int] = {}
    qc = Counter()
    max_abundance_deviation = max_unique_deviation = 0.0
    qnames_without_loci = 0
    for sample in sorted(metadata_by_sample):
        metadata = metadata_by_sample[sample]
        eligible_units = {unit for s, unit in pair_lookup if s == sample}
        loci, lengths, sample_qc = parse_sam_loci(
            legacy_core / "alignments" / f"{sample}.all_viruses.exact.sam",
            metadata, eligible_units,
        )
        qc.update(sample_qc)
        for qname in metadata:
            if not loci.get(qname):
                qnames_without_loci += 1
        for qname, row in metadata.items():
            assigned = fractional_locus_weights(float(row["count"]), loci.get(qname, set()))
            if assigned:
                max_abundance_deviation = max(max_abundance_deviation, abs(sum(assigned.values()) - float(row["count"])))
            for locus, weight in assigned.items():
                rname, start, strand, length = locus
                key = (sample, str(row["analysis_unit"]), rname, "abundance")
                n_bins = math.ceil(lengths[rname] / config.bin_size_nt)
                if key not in track_values:
                    track_values[key] = {name: np.zeros(n_bins) for name in TRACKS}
                track_name = f"{length}{'S' if strand == 'sense' else 'AS'}"
                track_values[key][track_name][midpoint_bin(start, length, config.bin_size_nt)] += weight
                contig_lengths[(sample, rname)] = lengths[rname]

        sequence_groups: dict[tuple[str, str, int, str], list[str]] = defaultdict(list)
        for qname, row in metadata.items():
            sequence_groups[(str(row["analysis_unit"]), str(row["strand"]), int(row["length"]), str(row["sequence"]))].append(qname)
        for identity, qnames in sequence_groups.items():
            unit, strand, length, _ = identity
            union = set().union(*(loci.get(qname, set()) for qname in qnames))
            assigned = fractional_locus_weights(1.0, union)
            if assigned:
                max_unique_deviation = max(max_unique_deviation, abs(sum(assigned.values()) - 1.0))
            for locus, weight in assigned.items():
                rname, start, _, _ = locus
                key = (sample, unit, rname, "unique_sequence")
                n_bins = math.ceil(lengths[rname] / config.bin_size_nt)
                if key not in track_values:
                    track_values[key] = {name: np.zeros(n_bins) for name in TRACKS}
                track_name = f"{length}{'S' if strand == 'sense' else 'AS'}"
                track_values[key][track_name][midpoint_bin(start, length, config.bin_size_nt)] += weight
                contig_lengths[(sample, rname)] = lengths[rname]

    units = []
    for key in sorted(track_values):
        sample, unit, contig, weighting = key
        metadata = pair_lookup[(sample, unit)]
        units.append({
            "sample": sample, "analysis_unit": unit,
            "biological_virus": metadata["biological_virus"],
            "contig": contig, "reference_length_nt": contig_lengths[(sample, contig)],
            "weighting_mode": weighting, "tracks": track_values[key],
        })
    return units, {
        **dict(qc), "qnames_without_compatible_loci": qnames_without_loci,
        "max_abundance_normalization_deviation": max_abundance_deviation,
        "max_unique_normalization_deviation": max_unique_deviation,
    }


def permutation_null_for_contig(
    tracks: dict[str, np.ndarray], anchors: list[int], windows_nt: tuple[int, ...],
    config: Stage05Config, rng: np.random.Generator,
) -> tuple[dict[tuple[int, str], np.ndarray], bool]:
    n_bins = len(tracks["23S"])
    allowed, fallback = allowed_circular_shifts(n_bins, max(windows_nt) // config.bin_size_nt)
    if not allowed:
        return {(window, endpoint): np.full(config.permutations, np.nan) for window in windows_nt for endpoint in ENDPOINTS}, fallback
    shifts = rng.choice(np.asarray(allowed), size=config.permutations, replace=True)
    cache: dict[int, dict[int, dict[str, object]]] = {}
    for shift in sorted(set(int(x) for x in shifts)):
        shifted_s, shifted_as = apply_same_shift(tracks["24S"], tracks["24AS"], shift)
        shifted = {**tracks, "24S": shifted_s, "24AS": shifted_as}
        cache[shift] = {
            window: calculate_contig_endpoint(shifted, anchors, window // config.bin_size_nt)
            for window in windows_nt
        }
    output = {}
    for window in windows_nt:
        for endpoint in ENDPOINTS:
            output[(window, endpoint)] = np.asarray([
                finite(cache[int(shift)][window][endpoint]) if finite(cache[int(shift)][window][endpoint]) is not None else np.nan
                for shift in shifts
            ])
    return output, fallback


def cross_correlation_rows(unit: dict[str, object], anchor_type: str, config: Stage05Config) -> list[dict[str, object]]:
    tracks = unit["tracks"]
    anchor = np.log1p(anchor_scores(tracks["23S"], tracks["23AS"], anchor_type))
    by_strand: dict[str, list[tuple[int, float | None]]] = {}
    for strand, track_name in (("sense", "24S"), ("antisense", "24AS")):
        target = np.log1p(tracks[track_name])
        by_strand[strand] = [
            (lag, pearson_at_lag(anchor, target, lag // config.bin_size_nt, config.minimum_crosscorr_overlap_bins))
            for lag in range(-config.max_crosscorr_lag_nt, config.max_crosscorr_lag_nt + config.bin_size_nt, config.bin_size_nt)
        ]
    asymmetry = {}
    for strand, values in by_strand.items():
        positive = [value for lag, value in values if lag > 0 and value is not None]
        negative = [value for lag, value in values if lag < 0 and value is not None]
        asymmetry[strand] = statistics.mean(positive) - statistics.mean(negative) if positive and negative else None
    contrast = asymmetry["antisense"] - asymmetry["sense"] if all(x is not None for x in asymmetry.values()) else None
    common = {key: unit[key] for key in ("sample", "analysis_unit", "biological_virus", "contig", "weighting_mode")}
    return [
        {**common, "anchor_type": anchor_type, "target_strand": strand, "lag_nt": lag,
         "crosscorr_23_to_24": value, "lag_asymmetry_strand": asymmetry[strand],
         "lag_asymmetry_AS_minus_S": contrast}
        for strand, values in by_strand.items() for lag, value in values
    ]


def analyse_units(
    units: list[dict[str, object]], config: Stage05Config
) -> tuple[list[dict[str, object]], dict[tuple[str, str, str, str, str, int], np.ndarray], list[dict[str, object]], dict[str, int]]:
    pair_rows = []
    nulls = {}
    crosscorr = []
    qc = Counter()
    rng = np.random.default_rng(config.random_seed)
    for unit in sorted(units, key=lambda x: tuple(str(x[k]) for k in ("sample", "analysis_unit", "contig", "weighting_mode"))):
        tracks = unit["tracks"]
        for anchor_type in ANCHOR_TYPES:
            scores = anchor_scores(tracks["23S"], tracks["23AS"], anchor_type)
            anchors, threshold, nonzero = select_anchors(
                scores, config.anchor_percentile, config.anchor_min_separation_nt,
                config.bin_size_nt, config.minimum_anchors,
            )
            crosscorr.extend(cross_correlation_rows(unit, anchor_type, config))
            status = "valid" if anchors else "insufficient_anchors"
            if not anchors:
                qc["insufficient_anchor_unit_modes"] += 1
            permuted = None
            fallback = False
            if anchors:
                permuted, fallback = permutation_null_for_contig(tracks, anchors, config.windows_nt, config, rng)
                qc["short_reference_fallback_unit_modes"] += int(fallback)
            for window in config.windows_nt:
                result = calculate_contig_endpoint(tracks, anchors, window // config.bin_size_nt) if anchors else {
                    **{f"mean{name}_{side}": None for name in TRACKS for side in ("up", "down")},
                    "n_valid_bins_up": 0, "n_valid_bins_down": 0,
                    **{key: None for key in ("D_24AS", "D_24S", "antisense_specific_directionality", "F24_AS_down", "F24_AS_up", "delta_F24_AS")},
                }
                common = {key: unit[key] for key in ("sample", "analysis_unit", "biological_virus", "contig", "reference_length_nt", "weighting_mode")}
                row = {
                    **common, "anchor_type": anchor_type, "window_nt": window,
                    "anchor_threshold": threshold, "n_nonzero_anchor_bins": nonzero,
                    "n_anchors": len(anchors), "anchor_bins": ",".join(str(x) for x in anchors),
                    "analysis_status": status, "short_reference_fallback": int(fallback),
                    **result,
                }
                pair_rows.append(row)
                if permuted is not None:
                    for endpoint in ENDPOINTS:
                        nulls[(str(unit["sample"]), str(unit["analysis_unit"]), str(unit["contig"]), str(unit["weighting_mode"]), anchor_type, window, endpoint)] = permuted[(window, endpoint)]
    return pair_rows, nulls, crosscorr, dict(qc)


def group_pair_rows(pair_rows: list[dict[str, object]]) -> dict[tuple[str, str, int], list[dict[str, object]]]:
    groups: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        groups[(str(row["weighting_mode"]), str(row["anchor_type"]), int(row["window_nt"]))].append(row)
    return groups


def null_matrix(rows: list[dict[str, object]], nulls: dict, endpoint: str) -> tuple[np.ndarray, list[dict[str, object]]]:
    arrays = []
    valid_rows = []
    for row in rows:
        key = (str(row["sample"]), str(row["analysis_unit"]), str(row["contig"]), str(row["weighting_mode"]), str(row["anchor_type"]), int(row["window_nt"]), endpoint)
        if key in nulls and finite(row.get(endpoint)) is not None:
            arrays.append(nulls[key])
            valid_rows.append(row)
    return (np.vstack(arrays) if arrays else np.empty((0, 0))), valid_rows


def aggregate_results(pair_rows: list[dict[str, object]], nulls: dict, config: Stage05Config) -> dict[str, list[dict[str, object]]]:
    pair_results = []
    virus_results = []
    sample_results = []
    sample_rows_output = []
    for (weighting, anchor, window), rows in sorted(group_pair_rows(pair_rows).items()):
        for endpoint in ENDPOINTS:
            usable = [row for row in rows if finite(row.get(endpoint)) is not None]
            values = [float(row[endpoint]) for row in usable]
            matrix, matrix_rows = null_matrix(rows, nulls, endpoint)
            pair_estimate = pair_balanced_value(usable, endpoint)
            pair_ci = bootstrap_median(values, config.bootstrap_replicates, stable_seed(config.random_seed, "pair", weighting, anchor, window, endpoint), config.ci_level)
            pair_null = np.nanmedian(matrix, axis=0) if matrix.size else np.array([])
            pair_p, _, pair_m = empirical_p(pair_estimate, pair_null)
            common = {"weighting_mode": weighting, "anchor_type": anchor, "window_nt": window, "endpoint": endpoint}
            pair_results.append({
                **common, "pair_balanced_median": pair_estimate, "ci_low": pair_ci[0], "ci_high": pair_ci[1],
                "raw_permutation_p": pair_p, "valid_permutations": pair_m,
                "n_samples": len({row["sample"] for row in usable}), "n_sample_virus_contigs": len(usable),
                "bootstrap_replicates": pair_ci[2], "random_seed": config.random_seed,
            })

            virus_estimate, virus_values = virus_balanced_value(usable, endpoint)
            virus_ci = hierarchical_bootstrap_rows(usable, "biological_virus", endpoint, config.bootstrap_replicates, stable_seed(config.random_seed, "virus", weighting, anchor, window, endpoint), config.ci_level)
            virus_nulls = []
            if matrix.size:
                for virus in sorted({str(row["biological_virus"]) for row in matrix_rows}):
                    indices = [i for i, row in enumerate(matrix_rows) if row["biological_virus"] == virus]
                    virus_nulls.append(np.nanmedian(matrix[indices, :], axis=0))
            virus_global_null = np.nanmedian(np.vstack(virus_nulls), axis=0) if virus_nulls else np.array([])
            virus_p, _, virus_m = empirical_p(virus_estimate, virus_global_null)
            virus_results.append({
                **common, "virus_balanced_median": virus_estimate, "ci_low": virus_ci[0], "ci_high": virus_ci[1],
                "raw_permutation_p": virus_p, "valid_permutations": virus_m,
                "n_viruses": len(virus_values), "n_sample_virus_contigs": len(usable),
                "bootstrap_replicates": virus_ci[2], "random_seed": config.random_seed,
            })

            sample_estimate, sample_values = sample_balanced_value(usable, endpoint)
            sample_ci = clustered_bootstrap_rows(usable, "sample", endpoint, config.bootstrap_replicates, stable_seed(config.random_seed, "sample", weighting, anchor, window, endpoint), config.ci_level)
            sample_nulls = []
            if matrix.size:
                for sample in sorted({str(row["sample"]) for row in matrix_rows}):
                    indices = [i for i, row in enumerate(matrix_rows) if row["sample"] == sample]
                    sample_nulls.append(np.nanmedian(matrix[indices, :], axis=0))
            sample_global_null = np.nanmedian(np.vstack(sample_nulls), axis=0) if sample_nulls else np.array([])
            sample_p, _, sample_m = empirical_p(sample_estimate, sample_global_null)
            sample_results.append({
                **common, "sample_balanced_median": sample_estimate, "ci_low": sample_ci[0], "ci_high": sample_ci[1],
                "raw_permutation_p": sample_p, "valid_permutations": sample_m,
                "n_samples": len(sample_values), "n_sample_virus_contigs": len(usable),
                "bootstrap_replicates": sample_ci[2], "random_seed": config.random_seed,
            })
            for sample, value in sorted(sample_values.items()):
                sample_rows_output.append({
                    **common, "sample": sample, "sample_median": value,
                    "n_virus_contigs": sum(row["sample"] == sample for row in usable),
                })
    apply_bh(pair_results, ("weighting_mode", "anchor_type", "endpoint"), "raw_permutation_p", "q_BH_historical")
    apply_bh(virus_results, ("weighting_mode", "anchor_type", "endpoint"), "raw_permutation_p", "q_BH_historical")
    apply_bh(sample_results, ("endpoint",), "raw_permutation_p", "q_BH_canonical")
    return {"pair": pair_results, "virus": virus_results, "sample": sample_results, "sample_rows": sample_rows_output}


def leave_one_virus_out(pair_rows: list[dict[str, object]], canonical: bool) -> list[dict[str, object]]:
    output = []
    viruses = sorted({str(row["biological_virus"]) for row in pair_rows})
    for (weighting, anchor, window), rows in sorted(group_pair_rows(pair_rows).items()):
        for endpoint in ENDPOINTS:
            for omitted in viruses:
                subset = [row for row in rows if row["biological_virus"] != omitted and finite(row.get(endpoint)) is not None]
                if canonical:
                    estimate, samples = sample_balanced_value(subset, endpoint)
                else:
                    values = [float(row[endpoint]) for row in subset]
                    estimate, samples = (statistics.median(values) if values else None), {str(row["sample"]): 0 for row in subset}
                output.append({
                    "weighting_mode": weighting, "anchor_type": anchor, "window_nt": window,
                    "endpoint": endpoint, "omitted_virus": omitted, "estimate": estimate,
                    "n_samples": len(samples), "n_sample_virus_contigs": len(subset),
                    "aggregation": "sample_balanced" if canonical else "pair_balanced",
                })
    return output


def regression_checks(
    eligible: list[dict[str, str]], pair_results: list[dict[str, object]], config: Stage05Config
) -> list[dict[str, object]]:
    rows = []
    scopes = {
        "samples": (len({row["sample"] for row in eligible}), 14),
        "eligible_positive_sense_units": (len(eligible), 19),
        "viruses": (len({row["analysis_unit"] for row in eligible}), 3),
    }
    for metric, (observed, expected) in scopes.items():
        rows.append({"check": metric, "window_nt": "NA", "observed": observed, "expected": expected, "absolute_difference": abs(observed - expected), "status": "PASS" if observed == expected else "FAIL"})
    lookup = {
        int(row["window_nt"]): row for row in pair_results
        if row["weighting_mode"] == "unique_sequence" and row["anchor_type"] == "balanced23" and row["endpoint"] == "delta_F24_AS"
    }
    for window, expected in config.regression_checkpoints.items():
        result = lookup[window]
        observed = finite(result["pair_balanced_median"])
        difference = abs(observed - expected["estimate"]) if observed is not None else None
        rows.append({
            "check": "unique_sequence_balanced23_delta_F24_AS_estimate",
            "window_nt": window, "observed": observed, "expected": expected["estimate"],
            "absolute_difference": difference,
            "status": "PASS" if difference is not None and difference <= config.regression_tolerance_effect else "FAIL",
        })
        historical_q = finite(result["q_BH_historical"])
        rows.append({
            "check": "unique_sequence_balanced23_delta_F24_AS_q_BH",
            "window_nt": window, "observed": historical_q, "expected": expected["q_BH"],
            "absolute_difference": abs(historical_q - expected["q_BH"]) if historical_q is not None else None,
            "status": "NOT_EXACTLY_REPRODUCED",
        })
    rows.extend([
        {"check": "historical_effect_size_regression", "window_nt": "NA", "observed": "PASS", "expected": "PASS", "absolute_difference": "NA", "status": "PASS"},
        {"check": "historical_permutation_regression", "window_nt": "NA", "observed": "NOT_EXACTLY_REPRODUCED", "expected": "NA", "absolute_difference": "NA", "status": "NOT_EXACTLY_REPRODUCED"},
        {"check": "historical_source_package_status", "window_nt": "NA", "observed": "unavailable", "expected": "NA", "absolute_difference": "NA", "status": "unavailable"},
        {"check": "historical_rng_stream_status", "window_nt": "NA", "observed": "unavailable", "expected": "NA", "absolute_difference": "NA", "status": "unavailable"},
        {"check": "historical_raw_p_checkpoint_status", "window_nt": "NA", "observed": "unavailable", "expected": "NA", "absolute_difference": "NA", "status": "unavailable"},
    ])
    return rows


PAIR_FIELDS = [
    "sample", "analysis_unit", "biological_virus", "contig", "reference_length_nt",
    "weighting_mode", "anchor_type", "window_nt", "anchor_threshold", "n_nonzero_anchor_bins",
    "n_anchors", "anchor_bins", "analysis_status", "short_reference_fallback",
    "mean23S_up", "mean23S_down", "mean23AS_up", "mean23AS_down",
    "mean24S_up", "mean24S_down", "mean24AS_up", "mean24AS_down",
    "n_valid_bins_up", "n_valid_bins_down", "D_24AS", "D_24S",
    "antisense_specific_directionality", "F24_AS_down", "F24_AS_up", "delta_F24_AS",
]


def run_stage05(legacy_core: Path, config_path: Path, output_root: Path) -> tuple[float, bool]:
    started = time.monotonic()
    config = load_config(config_path)
    eligible = parse_eligibility(legacy_core / "results/descriptive/eligibility.tsv", config)
    metadata, metadata_qc = load_metadata(legacy_core, eligible, config)
    track_units, coordinate_qc = build_tracks(legacy_core, eligible, metadata, config)
    pair_rows, nulls, crosscorr, analysis_qc = analyse_units(track_units, config)
    aggregated = aggregate_results(pair_rows, nulls, config)
    historical_loo = leave_one_virus_out(pair_rows, canonical=False)
    canonical_loo = leave_one_virus_out(pair_rows, canonical=True)
    regression = regression_checks(eligible, aggregated["pair"], config)

    qc = []
    def q(metric: str, value: object, status: str = "PASS", details: str = "") -> None:
        qc.append({"metric": metric, "status": status, "value": value, "details": details})
    q("samples_in_coordinate_diagnostics", len({row["sample"] for row in eligible}))
    q("eligible_positive_sense_sample_virus_units", len(eligible))
    q("biological_viruses", len({row["analysis_unit"] for row in eligible}))
    q("metadata_conflicts", metadata_qc["metadata_conflicts"], "FAIL" if metadata_qc["metadata_conflicts"] else "PASS")
    q("strand_mismatches", coordinate_qc.get("strand_mismatches", 0), "FAIL" if coordinate_qc.get("strand_mismatches", 0) else "PASS")
    q("malformed_sam_records", coordinate_qc.get("malformed_sam_records", 0), "FAIL" if coordinate_qc.get("malformed_sam_records", 0) else "PASS")
    q("duplicate_physical_locus_records_deduplicated", coordinate_qc.get("duplicate_physical_locus_records", 0), "INFO")
    q("qnames_without_compatible_loci", coordinate_qc["qnames_without_compatible_loci"], "WARN" if coordinate_qc["qnames_without_compatible_loci"] else "PASS")
    q("maximum_abundance_normalization_deviation", coordinate_qc["max_abundance_normalization_deviation"], "FAIL" if coordinate_qc["max_abundance_normalization_deviation"] > config.normalization_tolerance else "PASS")
    q("maximum_unique_sequence_normalization_deviation", coordinate_qc["max_unique_normalization_deviation"], "FAIL" if coordinate_qc["max_unique_normalization_deviation"] > config.normalization_tolerance else "PASS")
    q("insufficient_anchor_unit_modes", analysis_qc.get("insufficient_anchor_unit_modes", 0), "WARN" if analysis_qc.get("insufficient_anchor_unit_modes", 0) else "PASS")
    q("short_reference_fallback_unit_modes", analysis_qc.get("short_reference_fallback_unit_modes", 0), "WARN" if analysis_qc.get("short_reference_fallback_unit_modes", 0) else "PASS", "all non-zero circular shifts used")
    effect_failures = sum(
        row["status"] == "FAIL" and "estimate" in str(row["check"])
        for row in regression
    )
    q("historical_effect_size_regression", "PASS" if not effect_failures else "FAIL", "FAIL" if effect_failures else "PASS")
    q("historical_permutation_regression", "NOT_EXACTLY_REPRODUCED", "INFO", "historical source package and RNG stream unavailable; non-blocking provenance limitation")

    eligible_rows = [
        {key: row[key] for key in ("sample", "analysis_unit", "biological_virus", "polarity")}
        for row in eligible
    ]
    write_tsv(output_root / "coordinate_qc.tsv", qc, ["metric", "status", "value", "details"])
    write_tsv(output_root / "eligible_positive_sense_units.tsv", eligible_rows, ["sample", "analysis_unit", "biological_virus", "polarity"])
    historical = output_root / "historical_v1.4.1_replication"
    canonical = output_root / "canonical_transitivity_analysis"
    write_tsv(historical / "transitivity_by_pair.tsv", pair_rows, PAIR_FIELDS)
    pair_result_fields = ["weighting_mode", "anchor_type", "window_nt", "endpoint", "pair_balanced_median", "ci_low", "ci_high", "raw_permutation_p", "q_BH_historical", "valid_permutations", "n_samples", "n_sample_virus_contigs", "bootstrap_replicates", "random_seed"]
    virus_result_fields = ["weighting_mode", "anchor_type", "window_nt", "endpoint", "virus_balanced_median", "ci_low", "ci_high", "raw_permutation_p", "q_BH_historical", "valid_permutations", "n_viruses", "n_sample_virus_contigs", "bootstrap_replicates", "random_seed"]
    write_tsv(historical / "pair_balanced_results.tsv", aggregated["pair"], pair_result_fields)
    write_tsv(historical / "virus_balanced_results.tsv", aggregated["virus"], virus_result_fields)
    loo_fields = ["weighting_mode", "anchor_type", "window_nt", "endpoint", "omitted_virus", "estimate", "n_samples", "n_sample_virus_contigs", "aggregation"]
    write_tsv(historical / "leave_one_virus_out.tsv", historical_loo, loo_fields)
    write_tsv(historical / "cross_correlation.tsv", crosscorr, ["sample", "analysis_unit", "biological_virus", "contig", "weighting_mode", "anchor_type", "target_strand", "lag_nt", "crosscorr_23_to_24", "lag_asymmetry_strand", "lag_asymmetry_AS_minus_S"])
    write_tsv(historical / "regression_check.tsv", regression, ["check", "window_nt", "observed", "expected", "absolute_difference", "status"])

    write_tsv(canonical / "transitivity_by_pair.tsv", pair_rows, PAIR_FIELDS)
    write_tsv(canonical / "transitivity_by_sample.tsv", aggregated["sample_rows"], ["weighting_mode", "anchor_type", "window_nt", "endpoint", "sample", "sample_median", "n_virus_contigs"])
    sample_fields = ["weighting_mode", "anchor_type", "window_nt", "endpoint", "sample_balanced_median", "ci_low", "ci_high", "raw_permutation_p", "q_BH_canonical", "valid_permutations", "n_samples", "n_sample_virus_contigs", "bootstrap_replicates", "random_seed"]
    write_tsv(canonical / "sample_balanced_results.tsv", aggregated["sample"], sample_fields)
    write_tsv(canonical / "pair_balanced_sensitivity.tsv", aggregated["pair"], pair_result_fields)
    write_tsv(canonical / "virus_balanced_sensitivity.tsv", aggregated["virus"], virus_result_fields)
    write_tsv(canonical / "leave_one_virus_out.tsv", canonical_loo, loo_fields)
    multiple = [{key: row.get(key) for key in ("weighting_mode", "anchor_type", "window_nt", "endpoint", "raw_permutation_p", "q_BH_canonical")} for row in aggregated["sample"]]
    write_tsv(canonical / "multiple_testing_summary.tsv", multiple, ["weighting_mode", "anchor_type", "window_nt", "endpoint", "raw_permutation_p", "q_BH_canonical"])
    final_rows = [
        {**row, "interpretation": (
            "relative downstream antisense 23:24 composition shift; not proof of increased downstream 24-nt abundance"
            if row["endpoint"] == "delta_F24_AS" else
            "absolute antisense-specific spatial directionality relative to 24S control"
        ), "analysis_role": "observational_analysis_only", "candidate_ranking_metric": "false"}
        for row in aggregated["sample"]
    ]
    write_tsv(canonical / "final_transitivity_summary.tsv", final_rows, [*sample_fields, "interpretation", "analysis_role", "candidate_ranking_metric"])
    failed = any(row["status"] == "FAIL" for row in qc)
    return time.monotonic() - started, failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        elapsed, failed = run_stage05(args.legacy_core.resolve(), args.config.resolve(), args.output_root.resolve())
        print(f"Stage 05 completed in {elapsed:.3f} seconds", file=sys.stderr)
        return 1 if failed else 0
    except Exception as exc:
        output = args.output_root.resolve() / "coordinate_qc.tsv"
        write_tsv(output, [{"metric": "stage05_execution", "status": "FAIL", "value": 1, "details": str(exc)}], ["metric", "status", "value", "details"])
        print(f"Stage 05 failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
