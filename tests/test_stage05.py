import importlib.util
import dataclasses
import math
import sys
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).parents[1] / "workflow/scripts/stage05.py"
SPEC = importlib.util.spec_from_file_location("stage05", SCRIPT)
stage05 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = stage05
SPEC.loader.exec_module(stage05)


def config(permutations=50, bootstrap=50):
    return stage05.Stage05Config(
        ("+ssRNA",), (23, 24), 10, (100, 250, 500), 90.0, 50, 3,
        500, 10, permutations, bootstrap, 20260810, 0.95, 1e-12,
        1e-5, {},
    )


def test_alignment_midpoint_and_bin():
    assert stage05.alignment_midpoint_nt(10, 24) == 21.5
    assert stage05.alignment_midpoint_nt(10, 23) == 21
    assert stage05.midpoint_bin(10, 24) == 2


def test_abundance_split_deduplicates_loci_and_preserves_abundance():
    loci = [("c", 1), ("c", 1), ("c", 9)]
    weights = stage05.fractional_locus_weights(8, loci)
    assert weights == {("c", 1): 4, ("c", 9): 4}
    assert sum(weights.values()) == 8


def test_unique_sequence_across_k_loci_sums_exactly_one():
    weights = stage05.fractional_locus_weights(1, [("c", 1), ("c", 2), ("d", 1)])
    assert math.isclose(sum(weights.values()), 1.0)


def test_sam_duplicate_physical_locus_is_not_double_counted_and_mismatch_is_qc(tmp_path):
    sam = tmp_path / "x.sam"
    sam.write_text(
        "@SQ\tSN:V1|segment1\tLN:100\n"
        "Q1\t0\tV1|segment1\t11\t255\t23M\t*\t0\t0\t" + "A" * 23 + "\t*\n"
        "Q1\t0\tV1|segment1\t11\t255\t23M\t*\t0\t0\t" + "A" * 23 + "\t*\n"
        "Q2\t16\tV1|segment1\t21\t255\t23M\t*\t0\t0\t" + "C" * 23 + "\t*\n"
    )
    metadata = {
        "Q1": {"analysis_unit": "V1", "strand": "sense", "length": 23},
        "Q2": {"analysis_unit": "V1", "strand": "sense", "length": 23},
    }
    loci, lengths, qc = stage05.parse_sam_loci(sam, metadata, {"V1"})
    assert len(loci["Q1"]) == 1
    assert qc["duplicate_physical_locus_records"] == 1
    assert qc["strand_mismatches"] == 1
    assert lengths["V1|segment1"] == 100


def test_balanced_and_combined_anchor_formulas():
    sense = np.array([4.0, 1.0])
    antisense = np.array([9.0, 3.0])
    np.testing.assert_allclose(stage05.anchor_scores(sense, antisense, "balanced23"), [6, math.sqrt(3)])
    np.testing.assert_allclose(stage05.anchor_scores(sense, antisense, "combined23"), [13, 4])


def test_anchor_percentile_nonzero_tie_break_and_50nt_separation():
    scores = np.zeros(24)
    scores[[0, 6, 12]] = 10
    scores[[1, 2, 3, 4, 7, 8, 9]] = 1
    anchors, threshold, nonzero = stage05.select_anchors(scores, 90, 50, 10, 3)
    assert threshold == 10
    assert anchors == [0, 6, 12]
    assert nonzero == 10


def test_bins_exactly_50nt_apart_remain_within_historical_exclusion_radius():
    scores = np.zeros(20)
    scores[[0, 5, 10]] = 10
    scores[[1, 2, 3, 4, 6, 7, 8]] = 1
    anchors, _, _ = stage05.select_anchors(scores, 90, 50, 10, 3)
    assert anchors == []


def test_anchor_minimum_three_is_enforced():
    anchors, threshold, nonzero = stage05.select_anchors([0, 5, 0, 5], 90, 50, 10, 3)
    assert anchors == [] and threshold is None and nonzero == 2


def test_windows_exclude_anchor_use_direction_and_truncate_boundaries():
    track = np.array([0, 1, 2, 3, 4], dtype=float)
    upstream, n_up = stage05.pooled_window_mean(track, [0], 2, "upstream")
    downstream, n_down = stage05.pooled_window_mean(track, [0], 2, "downstream")
    assert upstream is None and n_up == 0
    assert downstream == 1.5 and n_down == 2


def test_overlapping_anchor_windows_repeat_genomic_bins():
    track = np.array([0, 1, 2, 3, 4], dtype=float)
    value, count = stage05.pooled_window_mean(track, [1, 2], 2, "downstream")
    assert value == 3 and count == 4  # values 2,3 and 3,4; bin 3 is repeated


def test_exact_endpoints_and_positive_composition_shift():
    means = {
        "mean24AS_down": 3.0, "mean24AS_up": 1.0,
        "mean24S_down": 2.0, "mean24S_up": 2.0,
        "mean23AS_down": 7.0, "mean23AS_up": 9.0,
    }
    result = stage05.endpoint_values(means)
    assert result["D_24AS"] == 0.5
    assert result["D_24S"] == 0
    assert result["antisense_specific_directionality"] == 0.5
    assert math.isclose(result["F24_AS_down"], 0.3)
    assert math.isclose(result["F24_AS_up"], 0.1)
    assert math.isclose(result["delta_F24_AS"], 0.2)


def test_endpoint_zero_denominators_are_na():
    means = {key: 0.0 for key in (
        "mean24AS_down", "mean24AS_up", "mean24S_down", "mean24S_up",
        "mean23AS_down", "mean23AS_up",
    )}
    result = stage05.endpoint_values(means)
    assert all(result[key] is None for key in ("D_24AS", "D_24S", "antisense_specific_directionality", "F24_AS_down", "F24_AS_up", "delta_F24_AS"))


def test_same_circular_shift_applies_to_both_24_tracks():
    sense = np.array([1, 2, 3, 4])
    antisense = np.array([10, 20, 30, 40])
    shifted_s, shifted_as = stage05.apply_same_shift(sense, antisense, 2)
    np.testing.assert_array_equal(shifted_s, [3, 4, 1, 2])
    np.testing.assert_array_equal(shifted_as, [30, 40, 10, 20])


def test_allowed_shift_exclusion_and_short_reference_fallback():
    allowed, fallback = stage05.allowed_circular_shifts(120, 50)
    assert not fallback and min(allowed) == 51 and max(allowed) == 69
    short, fallback = stage05.allowed_circular_shifts(80, 50)
    assert fallback and short == list(range(1, 80))


def test_empirical_one_sided_p_has_plus_one_correction():
    p, exceed, valid = stage05.empirical_p(2, [0, 2, 3, float("nan")])
    assert (p, exceed, valid) == (0.75, 2, 3)


def test_historical_pair_and_virus_medians():
    rows = [
        {"biological_virus": "V1", "value": 0},
        {"biological_virus": "V1", "value": 2},
        {"biological_virus": "V2", "value": 10},
    ]
    assert stage05.pair_balanced_value(rows, "value") == 2
    estimate, medians = stage05.virus_balanced_value(rows, "value")
    assert medians == {"V1": 1, "V2": 10}
    assert estimate == 5.5


def test_canonical_sample_median_of_virus_contigs():
    rows = [
        {"sample": "S1", "value": 0}, {"sample": "S1", "value": 2},
        {"sample": "S2", "value": 10},
    ]
    estimate, samples = stage05.sample_balanced_value(rows, "value")
    assert samples == {"S1": 1, "S2": 10}
    assert estimate == 5.5


def test_clustered_bootstrap_is_reproducible_and_keeps_sample_rows_together():
    rows = [
        {"sample": "S1", "value": 0}, {"sample": "S1", "value": 2},
        {"sample": "S2", "value": 10},
    ]
    first = stage05.clustered_bootstrap_rows(rows, "sample", "value", 100, 20260810, 0.95)
    second = stage05.clustered_bootstrap_rows(rows, "sample", "value", 100, 20260810, 0.95)
    assert first == second


def test_historical_three_window_and_canonical_twelve_test_bh_families():
    three = stage05.bh_adjust([0.01, 0.03, 0.2])
    assert np.allclose(three, [0.03, 0.045, 0.2])
    twelve = stage05.bh_adjust([0.01] * 12)
    assert len(twelve) == 12 and all(value == 0.01 for value in twelve)


def test_fixed_seed_permutation_reproducibility():
    cfg = config(permutations=20)
    tracks = {
        "23S": np.ones(120), "23AS": np.ones(120),
        "24S": np.arange(120, dtype=float), "24AS": np.arange(120, dtype=float)[::-1],
    }
    first, fallback1 = stage05.permutation_null_for_contig(tracks, [10, 40, 90], (100,), cfg, np.random.default_rng(cfg.random_seed))
    second, fallback2 = stage05.permutation_null_for_contig(tracks, [10, 40, 90], (100,), cfg, np.random.default_rng(cfg.random_seed))
    assert fallback1 == fallback2
    np.testing.assert_array_equal(first[(100, "delta_F24_AS")], second[(100, "delta_F24_AS")])


def test_historical_permutation_mismatch_is_recorded_but_effects_can_pass():
    cfg = dataclasses.replace(
        config(), regression_tolerance_effect=0.001,
        regression_checkpoints={100: {"estimate": 0.1, "q_BH": 0.2}},
    )
    pair = [{
        "weighting_mode": "unique_sequence", "anchor_type": "balanced23",
        "endpoint": "delta_F24_AS", "window_nt": 100,
        "pair_balanced_median": 0.1, "q_BH_historical": 0.7,
    }]
    eligible = [
        {"sample": f"S{i % 14}", "analysis_unit": ("V1", "V2", "V3")[i % 3]}
        for i in range(19)
    ]
    checks = stage05.regression_checks(eligible, pair, cfg)
    assert next(row for row in checks if row["check"] == "unique_sequence_balanced23_delta_F24_AS_estimate")["status"] == "PASS"
    assert next(row for row in checks if row["check"] == "unique_sequence_balanced23_delta_F24_AS_q_BH")["status"] == "NOT_EXACTLY_REPRODUCED"
