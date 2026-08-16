import importlib.util
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("stage09bc", REPO / "workflow/scripts/stage09bc.py")
stage09bc = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = stage09bc
SPEC.loader.exec_module(stage09bc)


@pytest.fixture(scope="module")
def transformed():
    source = REPO / "results/08_candidate_biophysics/candidate_biophysics.tsv"
    return stage09bc.transform(stage09bc.read_tsv(source))


def test_percentile_formula_ties_and_singleton():
    assert stage09bc.favourable_percentiles([1.0, 2.0, 2.0, 4.0]) == [0.125, 0.5, 0.5, 0.875]
    assert stage09bc.favourable_percentiles([9.0]) == [0.5]


def test_exact_candidate_row_preservation(transformed):
    assert len(transformed["layer2"]) == 1375
    assert len(transformed["layer3"]) == 1375
    assert Counter(int(row["candidate_length_nt"]) for row in transformed["layer2"]) == {23: 688, 24: 687}
    assert [row["candidate_id"] for row in transformed["layer2"]] == [row["candidate_id"] for row in transformed["layer3"]]


@pytest.mark.parametrize(
    "raw_column,percentile_column,layer",
    [
        ("asymmetry_ddg_4bp", "asymmetry_4bp_percentile", "layer2"),
        ("asymmetry_ddg_5bp", "asymmetry_5bp_percentile", "layer2"),
        ("guide_self_fold_mfe_kcal_mol", "guide_self_fold_percentile", "layer2"),
        ("target_whole_p_unpaired", "whole_site_accessibility_percentile", "layer3"),
        ("target_seed_g2_8_p_unpaired", "seed_accessibility_percentile", "layer3"),
    ],
)
def test_favourable_direction_is_higher(raw_column, percentile_column, layer, transformed):
    rows = [row for row in transformed[layer] if int(row["candidate_length_nt"]) == 23]
    low = min(rows, key=lambda row: float(row[raw_column]))
    high = max(rows, key=lambda row: float(row[raw_column]))
    assert float(low[percentile_column]) < float(high[percentile_column])


def test_normalization_is_separate_by_target_and_length(transformed):
    groups = defaultdict(list)
    for row in transformed["layer2"]:
        groups[(row["target_id"], int(row["candidate_length_nt"]))].append(row)
    assert {key[1]: len(rows) for key, rows in groups.items()} == {23: 688, 24: 687}
    for rows in groups.values():
        n = len(rows)
        percentiles = [float(row["asymmetry_4bp_percentile"]) for row in rows]
        assert min(percentiles) >= 0.5 / n
        assert max(percentiles) <= (n - 0.5) / n


def test_layer2_reference_score_uses_only_4bp_and_self_fold(transformed):
    for row in transformed["layer2"]:
        expected = 0.5 * float(row["asymmetry_4bp_percentile"]) + 0.5 * float(row["guide_self_fold_percentile"])
        assert float(row["layer2_reference_score"]) == pytest.approx(expected)
        assert float(row["asymmetry_4bp_5bp_percentile_difference"]) == pytest.approx(
            float(row["asymmetry_4bp_percentile"]) - float(row["asymmetry_5bp_percentile"])
        )


def test_layer3_reference_score_uses_canonical_whole_and_seed(transformed):
    for row in transformed["layer3"]:
        expected = 0.5 * float(row["whole_site_accessibility_percentile"]) + 0.5 * float(row["seed_accessibility_percentile"])
        assert float(row["layer3_reference_score"]) == pytest.approx(expected)


def test_w150_l100_is_canonical_and_alternatives_are_sensitivity_only(transformed):
    source = stage09bc.read_tsv(REPO / "results/08_candidate_biophysics/candidate_biophysics.tsv")
    indexes = [index for index, row in enumerate(source) if int(row["candidate_length_nt"]) == 23]
    expected_whole = stage09bc.favourable_percentiles([float(source[index]["target_whole_p_unpaired"]) for index in indexes])
    expected_seed = stage09bc.favourable_percentiles([float(source[index]["target_seed_g2_8_p_unpaired"]) for index in indexes])
    observed = [transformed["layer3"][index] for index in indexes]
    assert [float(row["whole_site_accessibility_percentile"]) for row in observed] == expected_whole
    assert [float(row["seed_accessibility_percentile"]) for row in observed] == expected_seed
    assert all("sensitivity_only" in row["status"] for row in transformed["layer3_sensitivity_23"])


def test_5bp_is_sensitivity_only_and_correlations_are_length_specific(transformed):
    assert all("sensitivity_only" in row["status"] for row in transformed["layer2_sensitivity_24"])
    correlations = transformed["layer2_correlations"]
    assert {(int(row["candidate_length_nt"]), row["comparison"]) for row in correlations} == {
        (23, "asymmetry_4bp_vs_self_fold"), (23, "asymmetry_4bp_vs_5bp"),
        (24, "asymmetry_4bp_vs_self_fold"), (24, "asymmetry_4bp_vs_5bp"),
    }
    assert all(row["spearman_rho"] is not None and math.isfinite(row["spearman_rho"]) for row in correlations)


def test_accessibility_sensitivity_correlations_are_complete(transformed):
    assert len(transformed["layer3_correlations"]) == 10
    assert Counter(int(row["candidate_length_nt"]) for row in transformed["layer3_correlations"]) == {23: 5, 24: 5}
    assert all(row["spearman_rho"] is not None and math.isfinite(row["spearman_rho"]) for row in transformed["layer3_correlations"])


def test_no_layer1_or_overall_score_fields(transformed):
    columns = set(transformed["layer2"][0]) | set(transformed["layer3"][0])
    assert not columns.intersection(stage09bc.FORBIDDEN_OUTPUT_COLUMNS)
    assert not any(column.startswith("guide_5p1_") for column in columns)
    assert "overall_score" not in columns and "overall_rank" not in columns
