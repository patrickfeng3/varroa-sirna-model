import importlib.util
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("stage10a", REPO / "workflow/scripts/stage10a.py")
stage10a = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = stage10a
SPEC.loader.exec_module(stage10a)


def identity(candidate_id, length, start=1):
    return {
        "target_id": "target_x",
        "transcript_id": "tx_x",
        "display_name": "Target X",
        "organism": "Synthetic",
        "candidate_id": candidate_id,
        "candidate_length_nt": str(length),
        "start_1based": str(start),
        "end_1based": str(start + length - 1),
        "target_sequence_dna": "A" * length,
        "target_sequence_rna": "A" * length,
        "antisense_guide_sequence_rna": "U" * length,
        "annotation_status": "unavailable",
        "start_region": "NA",
        "end_region": "NA",
        "overlap_regions": "NA",
        "crosses_annotation_boundary": "NA",
    }


def synthetic_inputs(specs):
    layer1, layer2, layer3 = [], [], []
    for index, (candidate_id, length, l1, l2, l3) in enumerate(specs, start=1):
        base = identity(candidate_id, length, index)
        layer1.append({**base, "layer1_accumulation_percentile": str(l1)})
        layer2.append({**base, "layer2_reference_score": str(l2)})
        layer3.append({**base, "layer3_reference_score": str(l3)})
    return layer1, layer2, layer3


@pytest.fixture(scope="module")
def current_result():
    layer1 = stage10a.read_tsv(REPO / "results/09_feature_layers/09A_layer1_accumulation/candidate_layer1.tsv")
    layer2 = stage10a.read_tsv(REPO / "results/09_feature_layers/09B_layer2_guide_competence/candidate_layer2.tsv")
    layer3 = stage10a.read_tsv(REPO / "results/09_feature_layers/09C_layer3_target_engagement/candidate_layer3.tsv")
    return (layer1, layer2, layer3, stage10a.transform(layer1, layer2, layer3, enforce_fixture=True))


def test_percentile_formula_ties_and_singleton():
    assert stage10a.favourable_percentiles([1.0, 2.0, 2.0, 4.0]) == [0.125, 0.5, 0.5, 0.875]
    assert stage10a.favourable_percentiles([7.0]) == [0.5]


def test_current_candidate_accounting_and_ids(current_result):
    layer1, layer2, layer3, result = current_result
    rows = result["candidates"]
    assert len(rows) == 1375
    assert Counter(int(row["candidate_length_nt"]) for row in rows) == {23: 688, 24: 687}
    assert [row["candidate_id"] for row in rows] == [row["candidate_id"] for row in layer1]
    assert {row["candidate_id"] for row in rows} == {row["candidate_id"] for row in layer2} == {
        row["candidate_id"] for row in layer3
    }


def test_identity_coordinates_and_sequences_preserved(current_result):
    layer1, _, _, result = current_result
    observed = {row["candidate_id"]: row for row in result["candidates"]}
    for source in layer1:
        row = observed[source["candidate_id"]]
        for column in stage10a.IDENTITY_COLUMNS:
            assert row[column] == source[column]


def test_layer1_is_exact_copy(current_result):
    layer1, _, _, result = current_result
    source = {row["candidate_id"]: float(row["layer1_accumulation_percentile"]) for row in layer1}
    assert all(row["stage10_layer1_percentile"] == source[row["candidate_id"]] for row in result["candidates"])


def test_layer2_and_layer3_percentiles_exact():
    inputs = synthetic_inputs([("a", 23, 0.2, 0.1, 0.4), ("b", 23, 0.4, 0.2, 0.3), ("c", 23, 0.6, 0.2, 0.2), ("d", 23, 0.8, 0.4, 0.1)])
    rows = stage10a.transform(*inputs)["candidates"]
    assert [row["stage10_layer2_percentile"] for row in rows] == [0.125, 0.5, 0.5, 0.875]
    assert [row["stage10_layer3_percentile"] for row in rows] == [0.875, 0.625, 0.375, 0.125]


def test_ranking_populations_are_separate_by_length():
    inputs = synthetic_inputs([("a23", 23, 0.1, 0.1, 0.1), ("b23", 23, 0.9, 0.9, 0.9), ("a24", 24, 0.2, 0.2, 0.2)])
    rows = {row["candidate_id"]: row for row in stage10a.transform(*inputs)["candidates"]}
    assert rows["a23"]["stage10_layer2_percentile"] == 0.25
    assert rows["b23"]["stage10_layer2_percentile"] == 0.75
    assert rows["a24"]["stage10_layer2_percentile"] == 0.5
    assert rows["a24"]["stage10_equal_layer_rank"] == 1.0


def test_equal_layer_score_rank_and_percentile_with_ties():
    inputs = synthetic_inputs([("a", 23, 0.5, 0.1, 0.1), ("b", 23, 0.5, 0.1, 0.1), ("c", 23, 0.9, 0.9, 0.9)])
    rows = {row["candidate_id"]: row for row in stage10a.transform(*inputs)["candidates"]}
    for row in rows.values():
        assert row["stage10_equal_layer_score"] == pytest.approx(
            (row["stage10_layer1_percentile"] + row["stage10_layer2_percentile"] + row["stage10_layer3_percentile"]) / 3
        )
    assert rows["a"]["stage10_equal_layer_rank"] == rows["b"]["stage10_equal_layer_rank"] == 2.5
    assert rows["a"]["stage10_equal_layer_percentile"] == rows["b"]["stage10_equal_layer_percentile"]
    assert rows["c"]["stage10_equal_layer_rank"] == 1.0


def test_pareto_dominance_and_fronts():
    vectors = [(1.0, 1.0, 1.0), (0.8, 0.9, 0.7), (0.9, 0.7, 0.8), (0.5, 0.5, 0.5)]
    assert stage10a.dominates(vectors[0], vectors[1])
    assert not stage10a.dominates(vectors[1], vectors[2])
    assert stage10a.pareto_fronts(vectors) == [1, 2, 2, 3]


def test_identical_vectors_do_not_dominate_and_share_front():
    vector = (0.8, 0.8, 0.8)
    assert not stage10a.dominates(vector, vector)
    assert stage10a.pareto_fronts([vector, vector, (0.2, 0.2, 0.2)]) == [1, 1, 2]


def test_every_current_candidate_has_one_positive_pareto_front(current_result):
    rows = current_result[3]["candidates"]
    assert len(rows) == len({row["candidate_id"] for row in rows})
    assert all(isinstance(row["stage10_pareto_front"], int) and row["stage10_pareto_front"] >= 1 for row in rows)


def test_minimum_layer_score_is_exact():
    inputs = synthetic_inputs([("a", 23, 0.2, 0.4, 0.8), ("b", 23, 0.8, 0.6, 0.2)])
    for row in stage10a.transform(*inputs)["candidates"]:
        assert row["stage10_minimum_layer_score"] == min(
            row["stage10_layer1_percentile"], row["stage10_layer2_percentile"], row["stage10_layer3_percentile"]
        )


def test_pareto_summary_is_exhaustive(current_result):
    result = current_result[3]
    expected = Counter((row["target_id"], int(row["candidate_length_nt"])) for row in result["candidates"])
    observed = defaultdict(int)
    fractions = defaultdict(float)
    for row in result["pareto_summary"]:
        key = (row["target_id"], int(row["candidate_length_nt"]))
        observed[key] += row["n_candidates"]
        fractions[key] += row["fraction_candidates"]
    assert observed == expected
    assert all(value == pytest.approx(1.0) for value in fractions.values())


def test_layer_correlations_are_three_per_stratum(current_result):
    result = current_result[3]
    grouped = Counter((row["target_id"], int(row["candidate_length_nt"])) for row in result["correlations"])
    assert set(grouped.values()) == {3}
    assert all(row["spearman_rho"] == "NA" or math.isfinite(row["spearman_rho"]) for row in result["correlations"])


def test_no_filter_gate_stage08_raw_or_stage11_outputs(current_result):
    layer1, _, _, result = current_result
    columns = set(result["candidates"][0])
    assert len(result["candidates"]) == len(layer1)
    assert not columns.intersection(stage10a.FORBIDDEN_STAGE08_RAW_COLUMNS)
    assert not columns.intersection(stage10a.FORBIDDEN_OUTPUT_COLUMNS)
    assert not any("region" in column.lower() for column in columns if column not in {"start_region", "end_region", "overlap_regions"})


def test_mismatched_candidate_identity_fails():
    inputs = list(synthetic_inputs([("a", 23, 0.5, 0.5, 0.5)]))
    inputs[1][0]["target_sequence_dna"] = "C" * 23
    with pytest.raises(ValueError, match="identity mismatch"):
        stage10a.transform(*inputs)
