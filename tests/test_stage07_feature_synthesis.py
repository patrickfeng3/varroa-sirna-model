"""Targeted tests for the post-hoc Stage 07 feature synthesis."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/stage07_feature_synthesis.py"
SPEC = importlib.util.spec_from_file_location("stage07_feature_synthesis", SCRIPT)
assert SPEC and SPEC.loader
SYNTHESIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SYNTHESIS)


def _positional_pair_rows(
    *, samples=("sample_1",), units=("virus_1",), lengths=(23,), strands=("antisense",)
):
    unique = {
        7: {"A": 0.18, "C": 0.19, "G": 0.21, "T": 0.42},
        10: {"A": 0.31, "C": 0.17, "G": 0.32, "T": 0.20},
        17: {"A": 0.29, "C": 0.23, "G": 0.22, "T": 0.26},
    }
    abundance = {
        7: {"A": 0.20, "C": 0.20, "G": 0.20, "T": 0.40},
        10: {"A": 0.35, "C": 0.15, "G": 0.35, "T": 0.15},
        17: {"A": 0.32, "C": 0.20, "G": 0.20, "T": 0.28},
    }
    rows = []
    for sample in samples:
        for unit in units:
            for length in lengths:
                for strand in strands:
                    for mode, frequencies in (
                        ("unique_sequence", unique), ("abundance", abundance)
                    ):
                        for position in (7, 10, 17):
                            for nucleotide in "ACGT":
                                rows.append({
                                    "sample": sample,
                                    "analysis_unit": unit,
                                    "biological_virus": unit,
                                    "polarity": "positive-sense",
                                    "length": str(length),
                                    "strand": strand,
                                    "weighting_mode": mode,
                                    "position_5p": str(position),
                                    "nucleotide": nucleotide,
                                    "observed_fraction": str(frequencies[position][nucleotide]),
                                    "expected_fraction": "0.25",
                                })
    return rows


def _grouped_pair_row(sample, unit, feature, length, strand, effect, delta):
    return {
        "sample": sample,
        "analysis_unit": unit,
        "biological_virus": unit,
        "polarity": "positive-sense",
        "length": length,
        "strand": strand,
        "feature_id": feature,
        "guide_position_5p": SYNTHESIS.GROUPED_FEATURES[feature]["position_5p"],
        "rna_bases": ",".join(SYNTHESIS.GROUPED_FEATURES[feature]["rna_bases"]),
        "tsv_bases": ",".join(SYNTHESIS.GROUPED_FEATURES[feature]["tsv_bases"]),
        "grouped_observed_fraction_unique": 0.5,
        "grouped_observed_fraction_abundance": 0.5,
        "grouped_expected_fraction": 0.5,
        "grouped_representation_enrichment_unique": effect,
        "grouped_representation_enrichment_abundance": effect,
        "grouped_representation_delta_unique": delta,
        "grouped_representation_delta_abundance": delta,
        "grouped_accumulation_ratio": effect,
        "grouped_accumulation_delta": delta,
    }


def _positional_summary_row(length, position, endpoint, effect):
    return {
        "length": str(length),
        "strand": "antisense",
        "endpoint": endpoint,
        "nucleotide": "A",
        "position_5p": str(position),
        "position_from_3p": str(length - position + 1),
        "sample_balanced_representation_enrichment": str(effect) if endpoint != "accumulation" else "NA",
        "sample_balanced_representation_delta_fraction": "0.1" if endpoint != "accumulation" else "NA",
        "sample_balanced_accumulation_ratio": str(effect) if endpoint == "accumulation" else "NA",
        "sample_balanced_accumulation_delta_fraction": "0.1" if endpoint == "accumulation" else "NA",
        "bootstrap_ci_low": "0.9",
        "bootstrap_ci_high": "1.2",
        "raw_p": "0.01",
        "bh_p": "0.02",
        "by_p": "0.03",
    }


def test_w7_group_is_a_plus_t():
    assert SYNTHESIS.GROUPED_FEATURES["W7"]["tsv_bases"] == ("A", "T")


def test_r10_group_is_a_plus_g():
    assert SYNTHESIS.GROUPED_FEATURES["R10"]["tsv_bases"] == ("A", "G")


def test_w17_group_is_a_plus_t():
    assert SYNTHESIS.GROUPED_FEATURES["W17"]["tsv_bases"] == ("A", "T")


def test_grouped_observed_and_expected_are_constituent_sums():
    rows, audits = SYNTHESIS.build_grouped_pair_rows(_positional_pair_rows())
    w7 = SYNTHESIS.find_row(rows, feature_id="W7")
    r10 = SYNTHESIS.find_row(rows, feature_id="R10")
    assert w7["grouped_observed_fraction_unique"] == pytest.approx(0.60)
    assert w7["grouped_observed_fraction_abundance"] == pytest.approx(0.60)
    assert r10["grouped_observed_fraction_unique"] == pytest.approx(0.63)
    assert w7["grouped_expected_fraction"] == pytest.approx(0.50)
    assert audits == {
        "max_expected_mode_difference": 0.0,
        "max_grouped_sum_difference": 0.0,
    }


def test_u_is_normalized_to_t():
    assert SYNTHESIS.normalize_rna_bases(("A", "U")) == ("A", "T")


def test_no_pseudocount_for_zero_denominator():
    assert SYNTHESIS.safe_ratio(1.0, 0.0) is None
    assert SYNTHESIS.safe_ratio(0.0, 0.0) is None
    assert SYNTHESIS.safe_ratio(0.0, 0.5) == 0.0


def test_exactly_18_antisense_family_tests():
    pair_rows = [
        _grouped_pair_row("sample_1", "virus_1", feature, length, strand, 1.1, 0.1)
        for feature in SYNTHESIS.GROUPED_FEATURES
        for length in SYNTHESIS.LENGTHS
        for strand in SYNTHESIS.STRANDS
    ]
    _, summary = SYNTHESIS.aggregate_grouped_features(pair_rows, bootstrap_replicates=20)
    assert len([row for row in summary if row["strand"] == "antisense"]) == 18
    assert all(
        row["antisense_family_bh_p"] is None
        for row in summary if row["strand"] == "sense"
    )


def test_pair_to_sample_to_dataset_median_ordering():
    pair_rows = [
        _grouped_pair_row("sample_1", "virus_1", "W7", 23, "antisense", 1.0, 0.1),
        _grouped_pair_row("sample_1", "virus_2", "W7", 23, "antisense", 3.0, 0.3),
        _grouped_pair_row("sample_2", "virus_3", "W7", 23, "antisense", 10.0, 0.5),
    ]
    samples, summary = SYNTHESIS.aggregate_grouped_features(pair_rows, bootstrap_replicates=20)
    sample_1 = SYNTHESIS.find_row(samples, sample="sample_1", feature_id="W7")
    accumulation = SYNTHESIS.find_row(summary, feature_id="W7", endpoint="accumulation")
    assert sample_1["grouped_accumulation_ratio"] == 2.0
    assert accumulation["sample_balanced_effect"] == 6.0
    assert accumulation["sample_balanced_delta"] == 0.35


def test_a3p3_coordinate_mapping_and_cross_length_support():
    rows = [
        _positional_summary_row(23, 21, "accumulation", 1.2),
        _positional_summary_row(24, 22, "accumulation", 1.3),
    ]
    cross = SYNTHESIS.build_cross_length(rows, "position_from_3p")
    row = SYNTHESIS.find_row(cross, relative_coordinate=3, nucleotide="A")
    assert SYNTHESIS.guide_position_5p_from_3p(23, 3) == 21
    assert SYNTHESIS.guide_position_5p_from_3p(24, 3) == 22
    assert row["length23_position_5p"] == "21"
    assert row["length24_position_5p"] == "22"
    assert row["supported_cross_length_match"] is True


def test_3p5_10_region_coordinate_mapping():
    assert SYNTHESIS.region_5p_from_3p(23, 5, 10) == (14, 19)
    assert SYNTHESIS.region_5p_from_3p(24, 5, 10) == (15, 20)


def test_regression_checkpoints_accept_expected_values():
    expected = {
        "W17": {23: (1.0202, 1.0648, 1.0404), 24: (1.0038, 1.0507, 1.0452)},
        "R10": {23: (1.0377, 1.0433, 0.9997), 24: (1.0100, 1.0173, 1.0167)},
        "W7": {23: (0.9857, 0.9938, 1.0108), 24: (0.9947, 0.9789, 0.9855)},
    }
    grouped = []
    endpoints = SYNTHESIS.ENDPOINTS
    for feature, length_values in expected.items():
        for length, effects in length_values.items():
            for endpoint, effect in zip(endpoints, effects):
                grouped.append({
                    "feature_id": feature, "length": length, "strand": "antisense",
                    "endpoint": endpoint, "sample_balanced_effect": effect,
                })
    positional = [
        _positional_summary_row(23, 21, "abundance_representation", 1.2373),
        _positional_summary_row(23, 21, "accumulation", 1.1358),
        _positional_summary_row(24, 22, "unique_representation", 1.1234),
        _positional_summary_row(24, 22, "abundance_representation", 1.3942),
        _positional_summary_row(24, 22, "accumulation", 1.2214),
    ]
    regional = [
        {"length": "23", "strand": "antisense", "region_3p": "GC_3p5-10",
         "endpoint": "accumulation", "sample_balanced_regional_gc6_delta": "-0.01623"},
        {"length": "24", "strand": "antisense", "region_3p": "GC_3p5-10",
         "endpoint": "accumulation", "sample_balanced_regional_gc6_delta": "-0.01949"},
    ]
    differences = SYNTHESIS.checkpoint_differences(grouped, positional, regional)
    assert max(differences.values()) <= 1e-12


@pytest.mark.parametrize(
    ("arguments", "expected"),
    [
        ({"external_prior": False, "compact_feature": True, "broad_context": False,
          "accumulation_corrected_both": True, "representation_positive_both": False,
          "representation_corrected_both": False, "predicted_direction_reproduced_both": True,
          "opposite_corrected_any": False}, "CARRY_FORWARD_HIGH"),
        ({"external_prior": True, "compact_feature": True, "broad_context": False,
          "accumulation_corrected_both": False, "representation_positive_both": True,
          "representation_corrected_both": True, "predicted_direction_reproduced_both": True,
          "opposite_corrected_any": False}, "CARRY_FORWARD_SUPPORTIVE"),
        ({"external_prior": False, "compact_feature": False, "broad_context": True,
          "accumulation_corrected_both": False, "representation_positive_both": True,
          "representation_corrected_both": True, "predicted_direction_reproduced_both": True,
          "opposite_corrected_any": False}, "CONTEXT_ONLY"),
        ({"external_prior": True, "compact_feature": True, "broad_context": False,
          "accumulation_corrected_both": False, "representation_positive_both": False,
          "representation_corrected_both": False, "predicted_direction_reproduced_both": False,
          "opposite_corrected_any": True}, "NOT_DEFAULT"),
    ],
)
def test_evidence_classification_is_rule_based(arguments, expected):
    assert SYNTHESIS.classify_evidence(**arguments) == expected


def test_evidence_schema_has_no_design_score_or_weight():
    forbidden = {"score", "weight", "bonus", "penalty", "rank"}
    assert forbidden.isdisjoint(SYNTHESIS.EVIDENCE_FIELDS)

