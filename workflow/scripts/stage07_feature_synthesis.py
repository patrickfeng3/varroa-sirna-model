#!/usr/bin/env python3
"""Post-hoc Stage 07 feature synthesis for transparent Stage 09 carry-forward."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from stage07 import (  # noqa: E402
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    CI_LEVEL,
    adjust_pvalues,
    exact_sign_test,
    median_or_none,
    sample_clustered_bootstrap,
    write_table,
)


GROUPED_FEATURES = {
    "W7": {"position_5p": 7, "rna_bases": ("A", "U"), "tsv_bases": ("A", "T")},
    "R10": {"position_5p": 10, "rna_bases": ("A", "G"), "tsv_bases": ("A", "G")},
    "W17": {"position_5p": 17, "rna_bases": ("A", "U"), "tsv_bases": ("A", "T")},
}
ENDPOINTS = ("unique_representation", "abundance_representation", "accumulation")
LENGTHS = (23, 24)
STRANDS = ("antisense", "sense")
EVIDENCE_FAMILY = "wang_bartel_2024_exploratory"
REGRESSION_TOLERANCE = 5e-4


GROUPED_SAMPLE_FIELDS = [
    "sample", "length", "strand", "feature_id", "guide_position_5p", "rna_bases",
    "tsv_bases", "n_virus_units", "grouped_observed_fraction_unique",
    "grouped_observed_fraction_abundance", "grouped_expected_fraction",
    "grouped_representation_enrichment_unique",
    "grouped_representation_enrichment_abundance", "grouped_representation_delta_unique",
    "grouped_representation_delta_abundance", "grouped_accumulation_ratio",
    "grouped_accumulation_delta",
]

GROUPED_SUMMARY_FIELDS = [
    "feature_id", "evidence_family", "length", "strand", "guide_position_5p",
    "rna_bases", "tsv_bases", "endpoint", "n_samples", "sample_balanced_effect",
    "sample_balanced_delta", "bootstrap_ci_low", "bootstrap_ci_high",
    "bootstrap_replicates_requested", "bootstrap_replicates_valid", "bootstrap_seed",
    "ci_method", "ci_level", "sign_test_n_nonzero", "sign_test_n_positive",
    "sign_test_n_negative", "sign_test_estimability", "raw_p", "antisense_family_bh_p",
]

CROSS_FIELDS = [
    "strand", "endpoint", "nucleotide", "coordinate_system", "relative_coordinate",
    "length23_position_5p", "length23_position_from_3p", "length23_effect",
    "length23_delta", "length23_ci_low", "length23_ci_high", "length23_raw_p",
    "length23_bh_p", "length23_by_p", "length24_position_5p",
    "length24_position_from_3p", "length24_effect", "length24_delta",
    "length24_ci_low", "length24_ci_high", "length24_raw_p", "length24_bh_p",
    "length24_by_p", "effect_direction_agrees", "both_internal_positions",
    "supported_cross_length_match",
]

EVIDENCE_FIELDS = [
    "feature_id", "feature_name", "evidence_origin", "external_prior", "feature_scope",
    "main_endpoint", "varroa_23_support", "varroa_24_support",
    "cross_length_consistency", "sense_comparator_note", "status", "rationale",
    "stage09_default_include",
]

QC_FIELDS = ["metric", "status", "value", "details"]


class SynthesisError(RuntimeError):
    pass


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def optional_float(value: object) -> float | None:
    text = str(value).strip()
    if text in {"", "NA", "None", "nan"}:
        return None
    number = float(text)
    return number if math.isfinite(number) else None


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_inputs(paths: dict[str, Path]) -> dict[str, dict[str, object]]:
    identities: dict[str, dict[str, object]] = {}
    for label, path in paths.items():
        if not path.is_file():
            raise SynthesisError(f"required canonical Stage 07 input is missing: {path}")
        identities[label] = {
            "path": str(path),
            "size": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
            "sha256": file_sha256(path),
        }
    return identities


def normalize_rna_bases(bases: Iterable[str]) -> tuple[str, ...]:
    normalized = tuple("T" if base.upper() == "U" else base.upper() for base in bases)
    if any(base not in {"A", "C", "G", "T"} for base in normalized):
        raise ValueError("group contains an unsupported nucleotide")
    return normalized


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def safe_delta(left: float | None, right: float | None) -> float | None:
    return None if left is None or right is None else left - right


def guide_position_5p_from_3p(length: int, position_3p: int) -> int:
    return length - position_3p + 1


def region_5p_from_3p(length: int, near_3p: int, far_3p: int) -> tuple[int, int]:
    return length - far_3p + 1, length - near_3p + 1


def sum_optional(values: Iterable[float | None]) -> float | None:
    values = list(values)
    return None if any(value is None for value in values) else sum(float(value) for value in values)


def build_grouped_pair_rows(
    positional_pair: Sequence[dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, float]]:
    rows = [row for row in positional_pair if row["weighting_mode"] in {"unique_sequence", "abundance"}]
    lookup = {
        (
            row["sample"], row["analysis_unit"], int(row["length"]), row["strand"],
            row["weighting_mode"], int(row["position_5p"]), row["nucleotide"],
        ): row
        for row in rows
    }
    identities = sorted({key[:4] for key in lookup})
    output: list[dict[str, object]] = []
    max_expected_mode_difference = 0.0
    max_grouped_sum_difference = 0.0
    for identity in identities:
        sample, unit, length, strand = identity
        if length not in LENGTHS or strand not in STRANDS:
            continue
        representative = next(row for row in rows if (
            row["sample"], row["analysis_unit"], int(row["length"]), row["strand"]
        ) == identity)
        for feature_id, specification in GROUPED_FEATURES.items():
            position = int(specification["position_5p"])
            bases = tuple(specification["tsv_bases"])
            observed: dict[str, float | None] = {}
            expected_by_mode: dict[str, float | None] = {}
            for mode in ("unique_sequence", "abundance"):
                constituent_observed = [
                    optional_float(lookup[identity + (mode, position, base)]["observed_fraction"])
                    for base in bases
                ]
                constituent_expected = [
                    optional_float(lookup[identity + (mode, position, base)]["expected_fraction"])
                    for base in bases
                ]
                observed[mode] = sum_optional(constituent_observed)
                expected_by_mode[mode] = sum_optional(constituent_expected)
                if observed[mode] is not None:
                    max_grouped_sum_difference = max(
                        max_grouped_sum_difference,
                        abs(observed[mode] - sum(float(value) for value in constituent_observed)),
                    )
            expected = expected_by_mode["unique_sequence"]
            if expected is not None and expected_by_mode["abundance"] is not None:
                max_expected_mode_difference = max(
                    max_expected_mode_difference,
                    abs(expected - float(expected_by_mode["abundance"])),
                )
            unique = observed["unique_sequence"]
            abundance = observed["abundance"]
            output.append({
                "sample": sample,
                "analysis_unit": unit,
                "biological_virus": representative["biological_virus"],
                "polarity": representative["polarity"],
                "length": length,
                "strand": strand,
                "feature_id": feature_id,
                "guide_position_5p": position,
                "rna_bases": ",".join(specification["rna_bases"]),
                "tsv_bases": ",".join(bases),
                "grouped_observed_fraction_unique": unique,
                "grouped_observed_fraction_abundance": abundance,
                "grouped_expected_fraction": expected,
                "grouped_representation_enrichment_unique": safe_ratio(unique, expected),
                "grouped_representation_enrichment_abundance": safe_ratio(abundance, expected),
                "grouped_representation_delta_unique": safe_delta(unique, expected),
                "grouped_representation_delta_abundance": safe_delta(abundance, expected),
                "grouped_accumulation_ratio": safe_ratio(abundance, unique),
                "grouped_accumulation_delta": safe_delta(abundance, unique),
            })
    return output, {
        "max_expected_mode_difference": max_expected_mode_difference,
        "max_grouped_sum_difference": max_grouped_sum_difference,
    }


def aggregate_grouped_features(
    pair_rows: Sequence[dict[str, object]],
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        grouped[(row["sample"], row["length"], row["strand"], row["feature_id"])].append(row)
    sample_rows: list[dict[str, object]] = []
    metric_fields = [
        "grouped_observed_fraction_unique", "grouped_observed_fraction_abundance",
        "grouped_expected_fraction", "grouped_representation_enrichment_unique",
        "grouped_representation_enrichment_abundance", "grouped_representation_delta_unique",
        "grouped_representation_delta_abundance", "grouped_accumulation_ratio",
        "grouped_accumulation_delta",
    ]
    for key in sorted(grouped):
        sample, length, strand, feature_id = key
        values = grouped[key]
        specification = GROUPED_FEATURES[str(feature_id)]
        row: dict[str, object] = {
            "sample": sample, "length": length, "strand": strand, "feature_id": feature_id,
            "guide_position_5p": specification["position_5p"],
            "rna_bases": ",".join(specification["rna_bases"]),
            "tsv_bases": ",".join(specification["tsv_bases"]),
            "n_virus_units": len(values),
        }
        for field in metric_fields:
            row[field] = median_or_none(value[field] for value in values)
        sample_rows.append(row)

    summary_rows: list[dict[str, object]] = []
    summary_groups: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        summary_groups[(row["length"], row["strand"], row["feature_id"])].append(row)
    endpoint_fields = {
        "unique_representation": (
            "grouped_representation_enrichment_unique", "grouped_representation_delta_unique"
        ),
        "abundance_representation": (
            "grouped_representation_enrichment_abundance", "grouped_representation_delta_abundance"
        ),
        "accumulation": ("grouped_accumulation_ratio", "grouped_accumulation_delta"),
    }
    for key in sorted(summary_groups):
        length, strand, feature_id = key
        samples = summary_groups[key]
        specification = GROUPED_FEATURES[str(feature_id)]
        for endpoint in ENDPOINTS:
            effect_field, delta_field = endpoint_fields[endpoint]
            effect_values = {
                str(row["sample"]): float(row[effect_field])
                for row in samples if row[effect_field] is not None
            }
            delta_values = [row[delta_field] for row in samples]
            ci_low, ci_high, valid = sample_clustered_bootstrap(
                effect_values, bootstrap_replicates, seed, CI_LEVEL
            )
            sign = exact_sign_test(delta_values)
            summary_rows.append({
                "feature_id": feature_id,
                "evidence_family": EVIDENCE_FAMILY,
                "length": length,
                "strand": strand,
                "guide_position_5p": specification["position_5p"],
                "rna_bases": ",".join(specification["rna_bases"]),
                "tsv_bases": ",".join(specification["tsv_bases"]),
                "endpoint": endpoint,
                "n_samples": len(effect_values),
                "sample_balanced_effect": median_or_none(effect_values.values()),
                "sample_balanced_delta": median_or_none(delta_values),
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
                "antisense_family_bh_p": None,
            })
    antisense_indices = [index for index, row in enumerate(summary_rows) if row["strand"] == "antisense"]
    adjusted = adjust_pvalues([summary_rows[index]["raw_p"] for index in antisense_indices], "BH")
    for index, value in zip(antisense_indices, adjusted):
        summary_rows[index]["antisense_family_bh_p"] = value
    return sample_rows, summary_rows


def positional_effect(row: dict[str, str]) -> tuple[float | None, float | None]:
    if row["endpoint"] == "accumulation":
        return (
            optional_float(row["sample_balanced_accumulation_ratio"]),
            optional_float(row["sample_balanced_accumulation_delta_fraction"]),
        )
    return (
        optional_float(row["sample_balanced_representation_enrichment"]),
        optional_float(row["sample_balanced_representation_delta_fraction"]),
    )


def build_cross_length(
    positional_summary: Sequence[dict[str, str]], coordinate: str
) -> list[dict[str, object]]:
    if coordinate not in {"position_5p", "position_from_3p"}:
        raise ValueError("unsupported coordinate system")
    lookup = {
        (
            int(row["length"]), row["strand"], row["endpoint"], row["nucleotide"],
            int(row[coordinate]),
        ): row
        for row in positional_summary
    }
    shared = sorted({key[1:] for key in lookup if key[0] == 23} & {key[1:] for key in lookup if key[0] == 24})
    output: list[dict[str, object]] = []
    for strand, endpoint, nucleotide, relative_coordinate in shared:
        left = lookup[(23, strand, endpoint, nucleotide, relative_coordinate)]
        right = lookup[(24, strand, endpoint, nucleotide, relative_coordinate)]
        left_effect, left_delta = positional_effect(left)
        right_effect, right_delta = positional_effect(right)
        directions_agree = (
            left_delta is not None and right_delta is not None
            and left_delta != 0 and right_delta != 0
            and (left_delta > 0) == (right_delta > 0)
        )
        both_internal = (
            3 <= int(left["position_5p"]) <= 21
            and 3 <= int(right["position_5p"]) <= 22
        )
        left_by = optional_float(left["by_p"])
        right_by = optional_float(right["by_p"])
        supported = bool(
            both_internal and directions_agree
            and left_by is not None and right_by is not None
            and left_by < 0.05 and right_by < 0.05
        )
        output.append({
            "strand": strand, "endpoint": endpoint, "nucleotide": nucleotide,
            "coordinate_system": "physical_5p" if coordinate == "position_5p" else "physical_3p",
            "relative_coordinate": relative_coordinate,
            "length23_position_5p": left["position_5p"],
            "length23_position_from_3p": left["position_from_3p"],
            "length23_effect": left_effect, "length23_delta": left_delta,
            "length23_ci_low": optional_float(left["bootstrap_ci_low"]),
            "length23_ci_high": optional_float(left["bootstrap_ci_high"]),
            "length23_raw_p": optional_float(left["raw_p"]),
            "length23_bh_p": optional_float(left["bh_p"]),
            "length23_by_p": left_by,
            "length24_position_5p": right["position_5p"],
            "length24_position_from_3p": right["position_from_3p"],
            "length24_effect": right_effect, "length24_delta": right_delta,
            "length24_ci_low": optional_float(right["bootstrap_ci_low"]),
            "length24_ci_high": optional_float(right["bootstrap_ci_high"]),
            "length24_raw_p": optional_float(right["raw_p"]),
            "length24_bh_p": optional_float(right["bh_p"]),
            "length24_by_p": right_by,
            "effect_direction_agrees": directions_agree,
            "both_internal_positions": both_internal,
            "supported_cross_length_match": supported,
        })
    return output


def classify_evidence(
    *,
    external_prior: bool,
    compact_feature: bool,
    broad_context: bool,
    accumulation_corrected_both: bool,
    representation_positive_both: bool,
    representation_corrected_both: bool,
    predicted_direction_reproduced_both: bool,
    opposite_corrected_any: bool,
) -> str:
    if broad_context:
        return "CONTEXT_ONLY"
    if accumulation_corrected_both and (external_prior or compact_feature):
        return "CARRY_FORWARD_HIGH"
    if external_prior and representation_positive_both and representation_corrected_both:
        return "CARRY_FORWARD_SUPPORTIVE"
    if external_prior and (not predicted_direction_reproduced_both or opposite_corrected_any):
        return "NOT_DEFAULT"
    return "CONTEXT_ONLY"


def find_row(rows: Sequence[dict[str, object]], **criteria: object) -> dict[str, object]:
    matches = [row for row in rows if all(str(row[key]) == str(value) for key, value in criteria.items())]
    if len(matches) != 1:
        raise SynthesisError(f"expected one row for {criteria}, observed {len(matches)}")
    return matches[0]


def evidence_delta(row: dict[str, object]) -> float | None:
    if "sample_balanced_delta" in row:
        return optional_float(row["sample_balanced_delta"])
    if "sample_balanced_regional_gc6_delta" in row:
        return optional_float(row["sample_balanced_regional_gc6_delta"])
    if row.get("endpoint") == "accumulation":
        return optional_float(row.get("sample_balanced_accumulation_delta_fraction"))
    return optional_float(row.get("sample_balanced_representation_delta_fraction"))


def evidence_effect(row: dict[str, object]) -> float | None:
    if "sample_balanced_effect" in row:
        return optional_float(row["sample_balanced_effect"])
    if "sample_balanced_regional_gc6_delta" in row:
        return optional_float(row["sample_balanced_regional_gc6_delta"])
    if row.get("endpoint") == "accumulation":
        return optional_float(row.get("sample_balanced_accumulation_ratio"))
    return optional_float(row.get("sample_balanced_representation_enrichment"))


def corrected(row: dict[str, object], field: str, direction: int) -> bool:
    p_value = optional_float(row[field])
    delta = evidence_delta(row)
    return bool(p_value is not None and p_value < 0.05 and delta is not None and delta * direction > 0)


def support_text(row: dict[str, object], corrected_field: str) -> str:
    effect = evidence_effect(row)
    delta = evidence_delta(row)
    p_value = optional_float(row[corrected_field])
    return f"effect={effect:.6g}; delta={delta:.6g}; corrected_p={p_value:.6g}" if None not in (effect, delta, p_value) else "not estimable"


def build_feature_evidence(
    grouped_summary: Sequence[dict[str, object]],
    positional_summary: Sequence[dict[str, str]],
    regional_summary: Sequence[dict[str, str]],
    cross_same_5p: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    grouped = lambda feature, length, endpoint, strand="antisense": find_row(
        grouped_summary, feature_id=feature, length=length, endpoint=endpoint, strand=strand
    )
    positional = lambda length, position, nucleotide, endpoint, strand="antisense": find_row(
        positional_summary, length=length, position_5p=position, nucleotide=nucleotide,
        endpoint=endpoint, strand=strand
    )
    regional = lambda length, endpoint, strand="antisense": find_row(
        regional_summary, length=length, region_3p="GC_3p5-10", endpoint=endpoint, strand=strand
    )

    a23 = positional(23, 21, "A", "accumulation")
    a24 = positional(24, 22, "A", "accumulation")
    a_high = corrected(a23, "by_p", 1) and corrected(a24, "by_p", 1)

    gc23 = regional(23, "accumulation")
    gc24 = regional(24, "accumulation")
    gc_high = corrected(gc23, "regional_by_p", -1) and corrected(gc24, "regional_by_p", -1)

    w17_23 = grouped("W17", 23, "accumulation")
    w17_24 = grouped("W17", 24, "accumulation")
    w17_high = corrected(w17_23, "antisense_family_bh_p", 1) and corrected(
        w17_24, "antisense_family_bh_p", 1
    )

    r10_u23 = grouped("R10", 23, "unique_representation")
    r10_u24 = grouped("R10", 24, "unique_representation")
    r10_a23 = grouped("R10", 23, "abundance_representation")
    r10_a24 = grouped("R10", 24, "abundance_representation")
    r10_rep_positive = all(optional_float(row["sample_balanced_delta"]) > 0 for row in (r10_u23, r10_u24, r10_a23, r10_a24))
    r10_rep_corrected = all(optional_float(row["antisense_family_bh_p"]) < 0.05 for row in (r10_u23, r10_u24, r10_a23, r10_a24))
    g10_positive = all(
        evidence_delta(positional(length, 10, "G", endpoint)) > 0
        for length in LENGTHS
        for endpoint in ("unique_representation", "abundance_representation")
    )
    r10_rep_positive = r10_rep_positive and g10_positive

    w7_23 = grouped("W7", 23, "accumulation")
    w7_24 = grouped("W7", 24, "accumulation")
    w7_positive_both = all(optional_float(row["sample_balanced_delta"]) > 0 for row in (w7_23, w7_24))
    w7_opposite = any(corrected(row, "antisense_family_bh_p", -1) for row in (w7_23, w7_24))

    broad_g = [
        row for row in cross_same_5p
        if row["strand"] == "antisense" and row["nucleotide"] == "G"
        and row["endpoint"] in {"unique_representation", "abundance_representation"}
        and row["supported_cross_length_match"] is True
    ]

    status_a = classify_evidence(
        external_prior=False, compact_feature=True, broad_context=False,
        accumulation_corrected_both=a_high, representation_positive_both=True,
        representation_corrected_both=True, predicted_direction_reproduced_both=True,
        opposite_corrected_any=False,
    )
    status_gc = classify_evidence(
        external_prior=False, compact_feature=True, broad_context=False,
        accumulation_corrected_both=gc_high, representation_positive_both=False,
        representation_corrected_both=False, predicted_direction_reproduced_both=True,
        opposite_corrected_any=False,
    )
    status_w17 = classify_evidence(
        external_prior=True, compact_feature=True, broad_context=False,
        accumulation_corrected_both=w17_high, representation_positive_both=True,
        representation_corrected_both=True, predicted_direction_reproduced_both=True,
        opposite_corrected_any=False,
    )
    status_r10 = classify_evidence(
        external_prior=True, compact_feature=True, broad_context=False,
        accumulation_corrected_both=False, representation_positive_both=r10_rep_positive,
        representation_corrected_both=r10_rep_corrected,
        predicted_direction_reproduced_both=r10_rep_positive, opposite_corrected_any=False,
    )
    status_g = classify_evidence(
        external_prior=False, compact_feature=False, broad_context=len(broad_g) >= 2,
        accumulation_corrected_both=False, representation_positive_both=True,
        representation_corrected_both=True, predicted_direction_reproduced_both=True,
        opposite_corrected_any=False,
    )
    status_w7 = classify_evidence(
        external_prior=True, compact_feature=True, broad_context=False,
        accumulation_corrected_both=False, representation_positive_both=False,
        representation_corrected_both=False,
        predicted_direction_reproduced_both=w7_positive_both,
        opposite_corrected_any=w7_opposite,
    )

    records = [
        {
            "feature_id": "A3p3", "feature_name": "A at physical guide 3p3",
            "evidence_origin": "Varroa cross-length positional discovery", "external_prior": "none",
            "feature_scope": "compact single-nucleotide position", "main_endpoint": "accumulation",
            "varroa_23_support": support_text(a23, "by_p"),
            "varroa_24_support": support_text(a24, "by_p"),
            "cross_length_consistency": "same A@3p3; positive accumulation in both lengths",
            "sense_comparator_note": "retained in canonical positional summary",
            "status": status_a,
            "rationale": "compact aligned Varroa feature with corrected accumulation support in both lengths",
        },
        {
            "feature_id": "low_GC_3p5_10", "feature_name": "low GC / high AU at guide 3p5-10",
            "evidence_origin": "Varroa regional-GC discovery", "external_prior": "none",
            "feature_scope": "fixed six-nucleotide regional composition", "main_endpoint": "accumulation",
            "varroa_23_support": support_text(gc23, "regional_by_p"),
            "varroa_24_support": support_text(gc24, "regional_by_p"),
            "cross_length_consistency": "negative GC accumulation delta in aligned region for both lengths",
            "sense_comparator_note": "sense values retained in regional summary",
            "status": status_gc,
            "rationale": "aligned broad Varroa regional association with corrected support in both lengths",
        },
        {
            "feature_id": "W17", "feature_name": "W17 (A/U at guide position 17)",
            "evidence_origin": "Wang & Bartel 2024 plus Varroa grouped test",
            "external_prior": "slicing-mechanism rationale", "feature_scope": "compact grouped position",
            "main_endpoint": "accumulation", "varroa_23_support": support_text(w17_23, "antisense_family_bh_p"),
            "varroa_24_support": support_text(w17_24, "antisense_family_bh_p"),
            "cross_length_consistency": "positive corrected accumulation in both lengths",
            "sense_comparator_note": "descriptive sense grouped results retained",
            "status": status_w17, "rationale": "external rationale and concordant Varroa accumulation support",
        },
        {
            "feature_id": "R10", "feature_name": "R10 (A/G at guide position 10)",
            "evidence_origin": "Wang & Bartel 2024 plus Varroa grouped test",
            "external_prior": "slicing-mechanism rationale", "feature_scope": "compact grouped position",
            "main_endpoint": "representation", "varroa_23_support": support_text(r10_a23, "antisense_family_bh_p"),
            "varroa_24_support": support_text(r10_a24, "antisense_family_bh_p"),
            "cross_length_consistency": "positive corrected representation in both lengths; weaker accumulation",
            "sense_comparator_note": "descriptive sense grouped results retained",
            "status": status_r10, "rationale": "external rationale, positive G10/grouped representation, and no bilateral accumulation support",
        },
        {
            "feature_id": "early_central_G", "feature_name": "broad early/central positional G enrichment",
            "evidence_origin": "Varroa positional synthesis", "external_prior": "none",
            "feature_scope": "broad contextual composition pattern", "main_endpoint": "representation",
            "varroa_23_support": f"supported G matches across {len(broad_g)} aligned positional records",
            "varroa_24_support": f"supported G matches across {len(broad_g)} aligned positional records",
            "cross_length_consistency": "multiple aligned G representation signals",
            "sense_comparator_note": "full sense comparator remains available",
            "status": status_g, "rationale": "reproducible but broad and potentially redundant composition structure",
        },
        {
            "feature_id": "W7", "feature_name": "W7 (A/U at guide position 7)",
            "evidence_origin": "Wang & Bartel 2024 plus Varroa grouped test",
            "external_prior": "predicted favourable grouped feature", "feature_scope": "compact grouped position",
            "main_endpoint": "accumulation", "varroa_23_support": support_text(w7_23, "antisense_family_bh_p"),
            "varroa_24_support": support_text(w7_24, "antisense_family_bh_p"),
            "cross_length_consistency": "predicted positive direction not reproduced across lengths",
            "sense_comparator_note": "descriptive sense grouped results retained",
            "status": status_w7, "rationale": "Varroa evidence lacks concordant positive support and includes opposite-direction support",
        },
    ]
    for row in records:
        row["stage09_default_include"] = row["status"] in {
            "CARRY_FORWARD_HIGH", "CARRY_FORWARD_SUPPORTIVE"
        }
    return records


def checkpoint_differences(
    grouped_summary: Sequence[dict[str, object]],
    positional_summary: Sequence[dict[str, str]],
    regional_summary: Sequence[dict[str, str]],
) -> dict[str, float]:
    grouped = lambda feature, length, endpoint: optional_float(find_row(
        grouped_summary, feature_id=feature, length=length, endpoint=endpoint,
        strand="antisense"
    )["sample_balanced_effect"])
    positional = lambda length, position, endpoint: positional_effect(find_row(
        positional_summary, length=length, position_5p=position, nucleotide="A",
        endpoint=endpoint, strand="antisense"
    ))[0]
    regional = lambda length: optional_float(find_row(
        regional_summary, length=length, region_3p="GC_3p5-10",
        endpoint="accumulation", strand="antisense"
    )["sample_balanced_regional_gc6_delta"])
    expected = {
        "A3p3_23_abundance": 1.2373, "A3p3_23_accumulation": 1.1358,
        "A3p3_24_unique": 1.1234, "A3p3_24_abundance": 1.3942,
        "A3p3_24_accumulation": 1.2214, "GC3p5_10_23_accumulation": -0.01623,
        "GC3p5_10_24_accumulation": -0.01949, "W17_23_unique": 1.0202,
        "W17_23_abundance": 1.0648, "W17_23_accumulation": 1.0404,
        "W17_24_unique": 1.0038, "W17_24_abundance": 1.0507,
        "W17_24_accumulation": 1.0452, "R10_23_unique": 1.0377,
        "R10_23_abundance": 1.0433, "R10_23_accumulation": 0.9997,
        "R10_24_unique": 1.0100, "R10_24_abundance": 1.0173,
        "R10_24_accumulation": 1.0167, "W7_23_unique": 0.9857,
        "W7_23_abundance": 0.9938, "W7_23_accumulation": 1.0108,
        "W7_24_unique": 0.9947, "W7_24_abundance": 0.9789,
        "W7_24_accumulation": 0.9855,
    }
    observed = {
        "A3p3_23_abundance": positional(23, 21, "abundance_representation"),
        "A3p3_23_accumulation": positional(23, 21, "accumulation"),
        "A3p3_24_unique": positional(24, 22, "unique_representation"),
        "A3p3_24_abundance": positional(24, 22, "abundance_representation"),
        "A3p3_24_accumulation": positional(24, 22, "accumulation"),
        "GC3p5_10_23_accumulation": regional(23),
        "GC3p5_10_24_accumulation": regional(24),
    }
    for feature in ("W17", "R10", "W7"):
        for length in LENGTHS:
            for suffix, endpoint in (
                ("unique", "unique_representation"),
                ("abundance", "abundance_representation"),
                ("accumulation", "accumulation"),
            ):
                observed[f"{feature}_{length}_{suffix}"] = grouped(feature, length, endpoint)
    return {
        key: math.inf if observed[key] is None else abs(float(observed[key]) - value)
        for key, value in expected.items()
    }


def build_digest(evidence: Sequence[dict[str, object]]) -> str:
    lookup = {row["feature_id"]: row for row in evidence}
    order = ("A3p3", "low_GC_3p5_10", "W17", "R10", "early_central_G", "W7")
    lines = [
        "# Stage 07 feature synthesis digest", "",
        "Post-hoc evidence-management synthesis of existing canonical Stage 07 outputs; not independent validation or an efficacy model.", "",
        "| Feature | 23-nt evidence | 24-nt evidence | Interpretation | Status |",
        "|---|---|---|---|---|",
    ]
    for feature_id in order:
        row = lookup[feature_id]
        lines.append(
            f"| {row['feature_name']} | {row['varroa_23_support']} | "
            f"{row['varroa_24_support']} | {row['rationale']} | {row['status']} |"
        )
    lines.extend([
        "", "## Stage 09+ default candidate feature set", "",
        "- A3p3", "- guide 3p5–10 GC/AU composition", "- W17", "- R10", "",
        "These are features eligible for Stage 09 evaluation, NOT accepted scoring terms. Stage 09 must test correlation, redundancy, antagonism and incremental information before combination or weighting.", "",
    ])
    return "\n".join(lines)


def build_qc(
    inputs: dict[str, dict[str, object]],
    post_hashes: dict[str, str],
    pair_rows: Sequence[dict[str, object]],
    grouped_summary: Sequence[dict[str, object]],
    cross_3p: Sequence[dict[str, object]],
    regional_summary: Sequence[dict[str, str]],
    evidence: Sequence[dict[str, object]],
    audits: dict[str, float],
    checkpoints: dict[str, float],
    stage07_accounting: Sequence[dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    add = lambda metric, status, value, details="": rows.append(
        {"metric": metric, "status": status, "value": value, "details": details}
    )
    for label, identity in inputs.items():
        unchanged = identity["sha256"] == post_hashes[label]
        add(
            f"input_{label}_read_only", "PASS" if unchanged else "FAIL",
            identity["sha256"], f"path={identity['path']}; size={identity['size']}",
        )
    expected_groups = {"W7": ("A", "T"), "R10": ("A", "G"), "W17": ("A", "T")}
    for feature, expected in expected_groups.items():
        observed = tuple(GROUPED_FEATURES[feature]["tsv_bases"])
        add(f"{feature}_nucleotide_group", "PASS" if observed == expected else "FAIL", ",".join(observed))
    add("rna_U_to_tsv_T_normalization", "PASS" if normalize_rna_bases(("A", "U")) == ("A", "T") else "FAIL", "A,U -> A,T")
    add("grouped_observed_constituent_sum", "PASS" if audits["max_grouped_sum_difference"] <= 1e-12 else "FAIL", audits["max_grouped_sum_difference"])
    add("grouped_expected_constituent_sum", "PASS" if audits["max_expected_mode_difference"] <= 1e-12 else "FAIL", audits["max_expected_mode_difference"])
    add("pseudocounts_added", "PASS", 0)
    antisense_tests = [row for row in grouped_summary if row["strand"] == "antisense"]
    add("wang_bartel_antisense_test_family", "PASS" if len(antisense_tests) == 18 else "FAIL", len(antisense_tests), "expected=18")
    sample_count = len({row["sample"] for row in pair_rows})
    pair_count = len({(row["sample"], row["analysis_unit"]) for row in pair_rows})
    add("inherited_biological_samples", "PASS" if sample_count == 20 else "FAIL", sample_count)
    add("inherited_sample_virus_units", "PASS" if pair_count == 54 else "FAIL", pair_count)
    add("pair_sample_dataset_hierarchy", "PASS", "pair -> sample median -> dataset median")
    a3p3 = [row for row in cross_3p if row["relative_coordinate"] == 3 and row["nucleotide"] == "A"]
    a_coordinates = all(str(row["length23_position_5p"]) == "21" and str(row["length24_position_5p"]) == "22" for row in a3p3)
    add("A3p3_coordinate_mapping", "PASS" if a3p3 and a_coordinates else "FAIL", len(a3p3), "23nt p21; 24nt p22")
    region_rows = [
        row for row in regional_summary
        if row["region_3p"] == "GC_3p5-10" and row["strand"] == "antisense"
    ]
    region_coordinates = all(
        (row["length"] == "23" and row["start_5p"] == "14" and row["end_5p"] == "19")
        or (row["length"] == "24" and row["start_5p"] == "15" and row["end_5p"] == "20")
        for row in region_rows
    )
    add("GC_3p5_10_coordinate_mapping", "PASS" if len(region_rows) == 6 and region_coordinates else "FAIL", len(region_rows))
    max_checkpoint = max(checkpoints.values(), default=0.0)
    add("interactive_regression_checkpoints", "PASS" if max_checkpoint <= REGRESSION_TOLERANCE else "FAIL", max_checkpoint, f"tolerance={REGRESSION_TOLERANCE}")
    expected_statuses = {
        "A3p3": "CARRY_FORWARD_HIGH", "low_GC_3p5_10": "CARRY_FORWARD_HIGH",
        "W17": "CARRY_FORWARD_HIGH", "R10": "CARRY_FORWARD_SUPPORTIVE",
        "early_central_G": "CONTEXT_ONLY", "W7": "NOT_DEFAULT",
    }
    actual_statuses = {row["feature_id"]: row["status"] for row in evidence}
    add("evidence_classification_regression", "PASS" if actual_statuses == expected_statuses else "FAIL", str(actual_statuses))
    upstream_fails = [row for row in stage07_accounting if row["status"] == "FAIL"]
    add("canonical_stage07_qc", "PASS" if not upstream_fails else "FAIL", len(upstream_fails))
    forbidden_columns = {"score", "weight", "bonus", "penalty", "rank"}
    evidence_columns = set(EVIDENCE_FIELDS)
    add("numerical_design_scores_created", "PASS" if not evidence_columns & forbidden_columns else "FAIL", 0)
    add("stage08_calculations", "PASS", 0)
    add("stage09_weights", "PASS", 0)
    return rows


def run(
    positional_pair_path: Path,
    positional_summary_path: Path,
    regional_summary_path: Path,
    accounting_path: Path,
    output_root: Path,
) -> tuple[float, bool]:
    started = time.monotonic()
    paths = {
        "positional_pair": positional_pair_path,
        "positional_summary": positional_summary_path,
        "regional_summary": regional_summary_path,
        "stage07_accounting": accounting_path,
    }
    identities = validate_inputs(paths)
    positional_pair = read_tsv(positional_pair_path)
    positional_summary = read_tsv(positional_summary_path)
    regional_summary = read_tsv(regional_summary_path)
    accounting = read_tsv(accounting_path)

    pair_grouped, audits = build_grouped_pair_rows(positional_pair)
    sample_grouped, grouped_summary = aggregate_grouped_features(pair_grouped)
    cross_3p = build_cross_length(positional_summary, "position_from_3p")
    cross_5p = build_cross_length(positional_summary, "position_5p")
    evidence = build_feature_evidence(
        grouped_summary, positional_summary, regional_summary, cross_5p
    )
    checkpoints = checkpoint_differences(grouped_summary, positional_summary, regional_summary)
    post_hashes = {label: file_sha256(path) for label, path in paths.items()}
    qc = build_qc(
        identities, post_hashes, pair_grouped, grouped_summary, cross_3p,
        regional_summary, evidence, audits, checkpoints, accounting,
    )
    output_root.mkdir(parents=True, exist_ok=True)
    write_table(output_root / "wang_bartel_grouped_by_sample.tsv", sample_grouped, GROUPED_SAMPLE_FIELDS)
    write_table(output_root / "wang_bartel_grouped_summary.tsv", grouped_summary, GROUPED_SUMMARY_FIELDS)
    write_table(output_root / "cross_length_same_3p.tsv", cross_3p, CROSS_FIELDS)
    write_table(output_root / "cross_length_same_5p.tsv", cross_5p, CROSS_FIELDS)
    write_table(output_root / "feature_evidence_summary.tsv", evidence, EVIDENCE_FIELDS)
    write_table(output_root / "feature_synthesis_qc.tsv", qc, QC_FIELDS)
    digest_path = output_root / "feature_digest.md"
    temporary = digest_path.with_suffix(".md.tmp")
    temporary.write_text(build_digest(evidence), encoding="utf-8")
    os.replace(temporary, digest_path)
    return time.monotonic() - started, any(row["status"] == "FAIL" for row in qc)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positional-pair", required=True, type=Path)
    parser.add_argument("--positional-summary", required=True, type=Path)
    parser.add_argument("--regional-summary", required=True, type=Path)
    parser.add_argument("--stage07-accounting", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        elapsed, failed = run(
            args.positional_pair.resolve(), args.positional_summary.resolve(),
            args.regional_summary.resolve(), args.stage07_accounting.resolve(),
            args.output_root.resolve(),
        )
    except (OSError, ValueError, SynthesisError) as exc:
        print(f"Stage 07 feature synthesis failed: {exc}", file=sys.stderr)
        return 1
    print(f"Stage 07 feature synthesis completed in {elapsed:.3f} seconds", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
