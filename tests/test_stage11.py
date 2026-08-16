import hashlib
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("stage11", REPO / "workflow/scripts/export_stage11_web_data.py")
stage11 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = stage11
SPEC.loader.exec_module(stage11)


@pytest.fixture(scope="module")
def current_export():
    payload, qc = stage11.build_payload(
        REPO / "results/10_candidate_integration/candidate_stage10a.tsv",
        REPO / "resources/targets/target_manifest.tsv",
        "Vd_CHIBIN",
        REPO,
        source_git_commit="test-commit",
        enforce_current_fixture=True,
    )
    return payload, qc


def synthetic_payload(lengths=(3, 4)):
    sequence = "ACGTACGT"
    candidates = []
    for length in lengths:
        for start in range(1, len(sequence) - length + 2):
            candidates.append({
                "candidate_id": f"x_{length}_{start}",
                "candidate_length_nt": length,
                "target_start_1based": start,
                "target_end_1based": start + length - 1,
                "layer1": start / 10,
                "layer2": (start + 1) / 10,
                "layer3": (start + 2) / 10,
                "total": (start + 3) / 10,
            })
    return {
        "schema_version": stage11.SCHEMA_VERSION,
        "target": {
            "target_id": "x", "display_name": "X", "transcript_length_nt": len(sequence),
            "transcript_sequence": sequence,
            "transcript_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "annotations": [
                {"feature": "5UTR", "start_1based": 1, "end_1based": 2},
                {"feature": "CDS", "start_1based": 3, "end_1based": 6},
                {"feature": "3UTR", "start_1based": 7, "end_1based": 8},
            ],
        },
        "supported_guide_lengths": list(lengths),
        "metrics": dict(stage11.METRIC_MAP),
        "candidates": candidates,
        "provenance": {},
    }


def test_export_preserves_current_ids_counts_and_scores(current_export):
    payload, qc = current_export
    source = stage11.read_tsv(REPO / "results/10_candidate_integration/candidate_stage10a.tsv")
    exported = {row["candidate_id"]: row for row in payload["candidates"]}
    assert set(exported) == {row["candidate_id"] for row in source}
    assert len(exported) == 1375
    assert Counter(row["candidate_length_nt"] for row in exported.values()) == {23: 688, 24: 687}
    for row in source:
        candidate = exported[row["candidate_id"]]
        assert candidate["layer1"] == float(row["stage10_layer1_percentile"])
        assert candidate["layer2"] == float(row["stage10_layer2_percentile"])
        assert candidate["layer3"] == float(row["stage10_layer3_percentile"])
        assert candidate["total"] == float(row["stage10_equal_layer_score"])
    assert all(row["status"] == "PASS" for row in qc)


def test_transcript_sequence_hash_annotations_and_schema(current_export):
    payload, _ = current_export
    target = payload["target"]
    assert payload["schema_version"] == "stage11-web-v1"
    assert len(target["transcript_sequence"]) == target["transcript_length_nt"] == 710
    assert hashlib.sha256(target["transcript_sequence"].encode()).hexdigest() == target["transcript_sequence_sha256"]
    assert target["annotations"] == [
        {"feature": "5UTR", "start_1based": 1, "end_1based": 329},
        {"feature": "CDS", "start_1based": 330, "end_1based": 665},
        {"feature": "3UTR", "start_1based": 666, "end_1based": 710},
    ]


def test_export_is_minimal_and_deterministic(current_export):
    payload, _ = current_export
    assert all(set(row) == stage11.EXPORTED_CANDIDATE_FIELDS for row in payload["candidates"])
    assert "stage10_pareto_front" not in stage11.compact_json(payload)
    assert "stage10_minimum_layer_score" not in stage11.compact_json(payload)
    assert stage11.compact_json(payload) == stage11.compact_json(json.loads(stage11.compact_json(payload)))
    javascript = stage11.compact_javascript(payload)
    assert javascript.startswith("window.STAGE11_PAYLOAD=")
    assert json.loads(javascript.removeprefix("window.STAGE11_PAYLOAD=").removesuffix(";\n")) == payload


def test_valid_region_and_contained_window_counts():
    payload = synthetic_payload()
    regions = stage11.compute_regions(payload, 3, 5, "layer1")
    assert len(regions) == 8 - 5 + 1
    assert all(row["stage11_region_n_contained_windows"] == 5 - 3 + 1 for row in regions)


@pytest.mark.parametrize("metric", ["layer1", "layer2", "layer3", "total"])
def test_r_equals_l_reproduces_exact_candidate_metric(metric):
    payload = synthetic_payload()
    regions = stage11.compute_regions(payload, 3, 3, metric)
    candidates = [row for row in payload["candidates"] if row["candidate_length_nt"] == 3]
    assert [row["stage11_region_mean_score"] for row in regions] == [row[metric] for row in candidates]


def test_invalid_region_lengths_produce_no_regions():
    payload = synthetic_payload()
    assert stage11.compute_regions(payload, 4, 3, "total") == []
    assert stage11.compute_regions(payload, 3, 9, "total") == []


def test_exact_region_mean_and_no_extra_smoothing():
    payload = synthetic_payload()
    regions = stage11.compute_regions(payload, 3, 5, "layer1")
    assert regions[0]["stage11_region_mean_score"] == pytest.approx((0.1 + 0.2 + 0.3) / 3)
    assert regions[1]["stage11_region_mean_score"] == pytest.approx((0.2 + 0.3 + 0.4) / 3)


def test_guide_lengths_never_mix():
    payload = synthetic_payload()
    regions3 = stage11.compute_regions(payload, 3, 5, "total")
    regions4 = stage11.compute_regions(payload, 4, 5, "total")
    assert regions3[0]["stage11_region_n_contained_windows"] == 3
    assert regions4[0]["stage11_region_n_contained_windows"] == 2
    assert regions3[0]["stage11_region_mean_score"] != regions4[0]["stage11_region_mean_score"]


def test_feature_boundary_classification():
    annotations = synthetic_payload()["target"]["annotations"]
    assert stage11.feature_for_position(annotations, 2) == "5UTR"
    assert stage11.feature_for_position(annotations, 3) == "CDS"
    assert stage11.feature_for_position(annotations, 6) == "CDS"
    assert stage11.feature_for_position(annotations, 7) == "3UTR"


def test_boundary_crossing_region_keeps_start_feature():
    regions = stage11.compute_regions(synthetic_payload(), 3, 5, "layer1")
    first = regions[0]
    assert first["stage11_region_start_1based"] == 1
    assert first["stage11_region_end_1based"] == 5
    assert first["stage11_region_start_feature"] == "5UTR"


def test_transcript_end_invalid_starts_are_excluded():
    regions = stage11.compute_regions(synthetic_payload(), 3, 5, "layer1")
    assert regions[-1]["stage11_region_start_1based"] == 4
    assert regions[-1]["stage11_region_end_1based"] == 8
    assert all(row["stage11_region_end_1based"] <= 8 for row in regions)


def test_top_bottom_sorting_overlap_and_padding():
    regions = stage11.compute_regions(synthetic_payload(), 3, 3, "layer1")
    tables = stage11.top_bottom_by_feature(regions)
    assert tables["5UTR"]["top"][0]["stage11_region_start_1based"] == 2
    assert tables["5UTR"]["bottom"][0]["stage11_region_start_1based"] == 1
    assert tables["5UTR"]["top"][2:] == [None, None, None]
    assert tables["3UTR"]["top"] == [None] * 5
    assert tables["CDS"]["top"][0]["stage11_region_start_1based"] - tables["CDS"]["top"][1]["stage11_region_start_1based"] == 1


def test_score_ties_use_lower_start_only_for_display_order():
    regions = stage11.compute_regions(synthetic_payload(), 3, 3, "layer1")
    for row in regions:
        row["stage11_region_mean_score"] = 0.5
    tables = stage11.top_bottom_by_feature(regions)
    starts = [row["stage11_region_start_1based"] for row in tables["CDS"]["top"] if row is not None]
    assert starts == sorted(starts)


def test_url_state_validation_and_defaults():
    assert stage11.validate_url_state({"guide": "4", "region": "5", "metric": "layer2"}, [3, 4], 8) == {
        "guide": 4, "region": 5, "metric": "layer2"
    }
    assert stage11.validate_url_state({"guide": "99", "region": "0", "metric": "bad"}, [3, 4], 8) == {
        "guide": 3, "region": 8, "metric": "total"
    }


def test_workflow_has_no_upstream_or_model_dependency():
    rule = (REPO / "workflow/rules/stage11.smk").read_text()
    assert "stage10=\"results/10_candidate_integration/candidate_stage10a.tsv\"" in rule
    assert "stage10=" in rule.split("params:", 1)[1]
    forbidden = ["stage00", "stage01", "stage02", "stage03", "stage04", "stage05", "stage06", "stage07", "stage08", "stage09", "ViennaRNA", "RNAfold", "RNAplfold", "model"]
    assert not any(token in rule for token in forbidden)


def test_guide_length_control_uses_selectable_radios():
    html = (REPO / "web/stage11/index.html").read_text()
    javascript = (REPO / "web/stage11/app.js").read_text()
    assert 'type="radio" name="guide-length" value="23"' in html
    assert 'type="radio" name="guide-length" value="24"' in html
    assert "renderGuideChoices(payload.supported_guide_lengths, state.guide)" in javascript
    assert "selectedGuideLength()" in javascript


def test_direct_file_payload_avoids_fetch_requirement():
    html = (REPO / "web/stage11/index.html").read_text()
    javascript = (REPO / "web/stage11/app.js").read_text()
    rule = (REPO / "workflow/rules/stage11.smk").read_text()
    assert '<script src="data/Vd_CHIBIN_stage11.js" defer></script>' in html
    assert "if (window.STAGE11_PAYLOAD)" in javascript
    assert "--web-data-js" in rule
