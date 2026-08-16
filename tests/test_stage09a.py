import csv
import gzip
import importlib.util
import math
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("stage09a", REPO / "workflow/scripts/stage09a.py")
stage09a = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = stage09a
SPEC.loader.exec_module(stage09a)


def encoded(sequence: str) -> dict[str, float]:
    return stage09a.encode_predictors(stage09a.guide_predictors(sequence))


def test_guide_orientation_and_all_eight_predictors():
    target = "ACGTTAGCTAGGCTAACGTACGT"
    guide = stage09a.reverse_complement_rna(target)
    assert guide == "ACGUACGUUAGCCUAGCUAACGU"
    values = stage09a.guide_predictors(guide)
    assert values == {
        "guide_5p1_nt": "A",
        "guide_5p2_nt": "C",
        "guide_3p2_nt": "G",
        "guide_3p1_nt": "U",
        "guide_A3p3": 0,
        "guide_GC_3p5_10": 2 / 6,
        "guide_W17": 0,
        "guide_R10": 1,
    }


def test_gc_3p5_10_is_exact_six_base_slice():
    guide = "A" * 13 + "GCGCAA" + "A" * 4
    assert len(guide) == 23
    assert guide[-10:-4] == "GCGCAA"
    assert stage09a.guide_predictors(guide)["guide_GC_3p5_10"] == 4 / 6


def test_terminal_encoding_uses_a_reference():
    features = encoded("A" * 23)
    assert tuple(features) == stage09a.BASE_FEATURE_NAMES
    assert all(features[name] == 0 for name in features if name.startswith("guide_") and name[-1] in "CGU")


def test_background_orientation_and_record_boundaries():
    records = ["A" * 23, "C" * 22]
    opportunities = stage09a.supported_antisense_sequences(records, 23)
    assert opportunities == {"T" * 23}
    assert "G" * 23 not in opportunities


def test_abundance_rows_are_summed(tmp_path):
    root = tmp_path
    table = root / "tables/S1/S1.read_level_features.tsv.gz"
    table.parent.mkdir(parents=True)
    fields = ["sample", "mapping_mode", "virus", "virus_assignment", "strand", "sequence", "length", "count"]
    with gzip.open(table, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for count in (2, 7):
            writer.writerow({"sample": "S1", "mapping_mode": "exact", "virus": "V1", "virus_assignment": "assigned",
                             "strand": "antisense", "sequence": "A" * 23, "length": 23, "count": count})
    observed, _ = stage09a.aggregate_observed_antisense(root, {("S1", "V1"): {}})
    assert observed[("S1", "V1", 23)]["A" * 23] == 9


def make_tiny_legacy(root: Path) -> None:
    eligibility = root / "results/descriptive/eligibility.tsv"
    eligibility.parent.mkdir(parents=True)
    eligibility.write_text(
        "sample\tanalysis_unit\tbiological_virus\tprimary_eligible\nS1\tV1\tFamily1\tTRUE\n",
        encoding="utf-8",
    )
    table = root / "tables/S1/S1.read_level_features.tsv.gz"
    table.parent.mkdir(parents=True)
    fields = ["sample", "mapping_mode", "virus", "virus_assignment", "strand", "sequence", "length", "count"]
    supported = "T" * 23
    with gzip.open(table, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"sample": "S1", "mapping_mode": "exact", "virus": "V1", "virus_assignment": "assigned",
                         "strand": "antisense", "sequence": supported, "length": 23, "count": 4})
        writer.writerow({"sample": "S1", "mapping_mode": "exact", "virus": "V1", "virus_assignment": "assigned",
                         "strand": "antisense", "sequence": "C" * 23, "length": 23, "count": 9})
    fasta = root / "references/consensus/S1.V1.final.background_masked.fa"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">r\n" + "A" * 24 + "\n", encoding="utf-8")


def test_outside_background_observation_is_excluded(tmp_path):
    make_tiny_legacy(tmp_path)
    accounting = stage09a.reconstruct_training_universe(tmp_path)
    assert accounting["represented_23nt"] == 1
    assert accounting["supported_abundance"] == 4
    assert accounting["outside_background_species"] == 1
    assert accounting["outside_background_abundance"] == 9


def test_exact_sample_aware_weights():
    rows = [
        {"sample": "S1", "analysis_unit": "V1", "candidate_length_nt": 23},
        {"sample": "S1", "analysis_unit": "V1", "candidate_length_nt": 23},
        {"sample": "S1", "analysis_unit": "V2", "candidate_length_nt": 23},
        {"sample": "S2", "analysis_unit": "V1", "candidate_length_nt": 23},
        {"sample": "S2", "analysis_unit": "V1", "candidate_length_nt": 23},
    ]
    weights = stage09a.sample_aware_weights(rows)
    assert math.isclose(sum(weights) / len(weights), 1)
    assert math.isclose(sum(weights[:3]), sum(weights[3:]))
    assert math.isclose(sum(weights[:2]), weights[2])


def test_fixed_effect_centering_is_weighted_within_group():
    x, y = stage09a.weighted_within_group_center(
        [[0.0], [2.0], [10.0], [14.0]], [1.0, 5.0, 3.0, 11.0],
        [1.0, 3.0, 2.0, 2.0], ["g1", "g1", "g2", "g2"],
    )
    for indexes, weights in (([0, 1], [1, 3]), ([2, 3], [2, 2])):
        assert sum(weights[i] * x[indexes[i]][0] for i in range(2)) == pytest.approx(0)
        assert sum(weights[i] * y[indexes[i]] for i in range(2)) == pytest.approx(0)


def test_cv_leakage_guard():
    stage09a.assert_cv_partition([{"biological_virus": "A"}], [{"biological_virus": "B"}], "biological_virus")
    with pytest.raises(ValueError, match="leakage"):
        stage09a.assert_cv_partition([{"biological_virus": "A"}], [{"biological_virus": "A"}], "biological_virus")


def test_top10_abundance_metrics_exact():
    scores = list(range(10))
    abundance = [1] * 9 + [9]
    share, lift = stage09a.top10_abundance_metrics(scores, abundance)
    assert share == 0.5
    assert lift == 5.0


def test_23_24_percentiles_are_separate_and_ties_average():
    rows = [(23, 1.0), (23, 3.0), (24, 100.0), (24, 100.0)]
    output = {}
    for length in (23, 24):
        indexes = [i for i, row in enumerate(rows) if row[0] == length]
        values = [rows[i][1] for i in indexes]
        for i, value in zip(indexes, stage09a.favourable_percentiles(values)):
            output[i] = value
    assert [output[i] for i in range(4)] == [0.25, 0.75, 0.5, 0.5]


def test_candidate_row_preservation_and_no_stage08_leakage(tmp_path):
    source = tmp_path / "candidates.tsv"
    target = "A" * 23
    fields = ["target_id", "candidate_id", "candidate_length_nt", "target_sequence_rna", "antisense_guide_sequence_rna"]
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"target_id": "generic", "candidate_id": "c1", "candidate_length_nt": 23,
                         "target_sequence_rna": target, "antisense_guide_sequence_rna": "U" * 23})
    destination = tmp_path / "prepared.tsv"
    assert stage09a.prepare_candidates(source, destination) == 1
    with destination.open() as handle:
        prepared = list(csv.DictReader(handle, delimiter="\t"))
    assert len(prepared) == 1 and prepared[0]["candidate_id"] == "c1"
    with pytest.raises(ValueError, match="forbidden"):
        stage09a.assert_no_stage08_feature_leakage(["asymmetry_ddg_4bp"])


def test_rscript_is_resolved_from_active_stage09a_environment(tmp_path):
    rscript = tmp_path / "bin/Rscript"
    rscript.parent.mkdir()
    rscript.write_text("", encoding="utf-8")
    assert stage09a.resolve_rscript({"CONDA_PREFIX": str(tmp_path)}) == rscript
    with pytest.raises(RuntimeError, match="CONDA_PREFIX"):
        stage09a.resolve_rscript({})


def test_v020_has_two_length_models_and_no_superseded_model_machinery():
    model = (REPO / "workflow/scripts/stage09a_model.R").read_text(encoding="utf-8")
    lowered = model.lower()
    assert "glmnet" not in lowered
    assert "alpha_grid" not in lowered
    assert "l1_ratio" not in lowered
    assert "nested" not in lowered
    assert "fit_representation" not in lowered
    assert "hurdle_score" not in lowered
    assert "lengths <- c(23l, 24l)" in lowered
    assert '"primary_accumulation_fit_count"' in model
    assert "12L" in model


def test_frozen_accounting_regression():
    config = REPO / "config/paths.local.yaml"
    if not config.is_file():
        pytest.skip("machine-local frozen core is not configured")
    legacy_line = next(line for line in config.read_text().splitlines() if line.strip().startswith("legacy_core:"))
    legacy = Path(legacy_line.split(":", 1)[1].strip().strip("'\""))
    if not legacy.is_dir():
        pytest.skip("configured frozen core is unavailable")
    accounting = stage09a.reconstruct_training_universe(legacy)
    stage09a.validate_frozen_accounting(accounting)
    assert {key: float(accounting[key]) for key in stage09a.EXPECTED_ACCOUNTING} == {
        key: float(value) for key, value in stage09a.EXPECTED_ACCOUNTING.items()
    }
