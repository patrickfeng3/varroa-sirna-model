import importlib.util
import math
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "workflow/scripts/stage04.py"
SPEC = importlib.util.spec_from_file_location("stage04", SCRIPT)
stage04 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = stage04
SPEC.loader.exec_module(stage04)


CONFIG = stage04.Stage04Config(100, 20260814, "percentile", 0.95, 1e-12, 2, -2)


def test_sample_balanced_aggregation_is_virus_then_sample_median():
    rows = [
        {"sample": "S1", "analysis_unit": "V1", "group": "G", "value": 0},
        {"sample": "S1", "analysis_unit": "V2", "group": "G", "value": 2},
        {"sample": "S2", "analysis_unit": "V1", "group": "G", "value": 10},
    ]
    sample, across = stage04.aggregate_metric_rows(rows, ("group",), "value", CONFIG)
    assert next(x for x in sample if x["sample"] == "S1")["sample_median"] == 1
    assert across[0]["sample_balanced_median"] == 5.5


def test_sample_clustered_bootstrap_is_reproducible():
    values = {"S1": 1.0, "S2": 2.0, "S3": 8.0}
    assert stage04.clustered_bootstrap(values, 1000, 20260814, 0.95) == stage04.clustered_bootstrap(
        values, 1000, 20260814, 0.95
    )


def spectrum_row(log_ratio=2, unique_log=1, count=999, z=100):
    return {
        "sample": "S1", "analysis_unit": "V1", "biological_virus": "V1",
        "focal_length": "23", "focal_strand": "sense", "end": "5p",
        "signed_distance": "0", "official_duplex_count": str(count),
        "official_unique_reference_count": "7",
        "official_duplex_log_ratio": str(log_ratio),
        "official_duplex_wald_z": str(z),
        "official_unique_reference_log_ratio": str(unique_log),
        "official_unique_reference_wald_z": "50", "run_id": "R",
    }


def test_full_spectrum_primary_effect_uses_official_log_ratio_not_count_or_z():
    _, _, across = stage04.aggregate_full_spectrum([spectrum_row()], CONFIG)
    duplex = next(x for x in across if x["official_view"] == "duplex")
    assert duplex["sample_balanced_steprna_log_ratio"] == 2
    assert duplex["sample_balanced_steprna_wald_z_descriptive"] == 100


def test_terminal_coordinates_preserve_physical_orientation():
    assert stage04.terminal_bases("ACAAAAATG") == {
        "5p1": "A", "5p2": "C", "3p2": "T", "3p1": "G"
    }


def focal_fixture():
    return [
        {
            "focal_id": "F1", "sequence": "AC" + "A" * 19 + "TG",
            "sample": "S1", "analysis_unit": "V1", "biological_virus": "V1",
            "focal_length": "23", "focal_strand": "antisense",
            "focal_abundance": "10", "run_id": "R1",
        },
        {
            "focal_id": "F2", "sequence": "GA" + "C" * 19 + "AT",
            "sample": "S1", "analysis_unit": "V1", "biological_virus": "V1",
            "focal_length": "23", "focal_strand": "antisense",
            "focal_abundance": "1", "run_id": "R1",
        },
        {
            "focal_id": "F3", "sequence": "TT" + "G" * 19 + "CC",
            "sample": "S1", "analysis_unit": "V1", "biological_virus": "V1",
            "focal_length": "23", "focal_strand": "antisense",
            "focal_abundance": "100", "run_id": "R1",
        },
    ]


def stage02_fixture():
    expected, general = [], []
    for mode in stage04.WEIGHTING_MODES:
        for position in stage04.POSITIONS:
            for nucleotide in stage04.NUCLEOTIDES:
                common = {
                    "sample": "S1", "analysis_unit": "V1", "length": "23",
                    "strand_scope": "antisense", "weighting_mode": mode,
                    "terminal_position": position, "nucleotide": nucleotide,
                }
                expected.append({**common, "expected_fraction": "0.25"})
                general.append({**common, "enrichment_ratio": "2"})
    return expected, general


def calculate(joint_ids=("F1",), recovered_ids=("F1", "F2")):
    joint = [
        {
            "focal_id": focal_id, "run_id": "R1",
            "focal_abundance": next(x["focal_abundance"] for x in focal_fixture() if x["focal_id"] == focal_id),
        }
        for focal_id in joint_ids
    ]
    expected, general = stage02_fixture()
    return stage04.calculate_sequence_pair_rows(
        focal_fixture(), joint, {"R1": set(recovered_ids)}, expected, general
    )


def find_pair(rows, mode, position="5p1", nucleotide="A"):
    return next(
        row for row in rows
        if row["weighting_mode"] == mode and row["terminal_position"] == position
        and row["nucleotide"] == nucleotide
    )


def test_unique_sequence_weighting_gives_each_reference_one():
    rows, _ = calculate()
    row = find_pair(rows, "unique_sequence")
    assert row["joint_observed_fraction"] == 1
    assert row["recovered_observed_fraction"] == 0.5


def test_abundance_weighting_uses_focal_abundance_not_passenger_multiplicity():
    rows, _ = calculate()
    row = find_pair(rows, "abundance")
    assert row["joint_total_weight"] == 10
    assert row["recovered_total_weight"] == 11
    assert math.isclose(row["recovered_observed_fraction"], 10 / 11)


def test_exact_stage02_background_and_general_enrichment_are_reused():
    rows, qc = calculate()
    row = find_pair(rows, "unique_sequence")
    assert row["stage02_expected_fraction"] == 0.25
    assert row["E_joint_absolute"] == 4
    assert row["E_all"] == 2
    assert qc["missing_expected"] == qc["missing_general"] == 0


def test_joint_vs_all_and_recovered_contrasts():
    rows, _ = calculate()
    unique = find_pair(rows, "unique_sequence")
    abundance = find_pair(rows, "abundance")
    assert unique["joint_vs_all_log2_contrast"] == 1
    assert unique["joint_vs_recovered_log2_contrast"] == 1
    assert math.isclose(abundance["joint_vs_recovered_log2_contrast"], math.log2(1.1))


def test_zero_required_fraction_produces_na_without_pseudocount():
    rows, _ = calculate()
    row = find_pair(rows, "unique_sequence", nucleotide="C")
    assert row["joint_observed_fraction"] == 0
    assert row["joint_vs_all_log2_contrast"] is None
    assert row["joint_vs_recovered_log2_contrast"] is None


def test_empty_joint_subset_produces_na_not_zero():
    rows, qc = calculate(joint_ids=())
    row = find_pair(rows, "unique_sequence")
    assert row["joint_observed_fraction"] is None
    assert row["E_joint_absolute"] is None
    assert qc["empty_joint"] == 1


def test_joint_subset_must_be_passenger_recovered():
    with pytest.raises(stage04.Stage04Error, match="outside recovered subset"):
        calculate(joint_ids=("F1",), recovered_ids=("F2",))


def test_recovered_ids_must_exist_in_focal_manifest():
    with pytest.raises(stage04.Stage04Error, match="absent from focal manifest"):
        calculate(joint_ids=(), recovered_ids=("UNKNOWN",))


def complete_population_rows():
    recovery, joint = [], []
    for sample, unit, d23, d24 in (
        ("S1", "V1", 0.1, 0.2),
        ("S1", "V2", 0.2, 0.5),
        ("S2", "V1", 0.1, 0.9),
    ):
        for length, value in ((23, d23), (24, d24)):
            base = {
                "sample": sample, "analysis_unit": unit, "focal_length": str(length),
                "focal_strand": "antisense",
            }
            recovery.append({
                **base,
                "passenger_recovery_fraction_unique": str(value),
                "passenger_recovery_fraction_abundance": str(value),
            })
            joint.append({**base, **{metric: str(value) for metric in stage04.JOINT_METRICS}})
    return recovery, joint


def test_paired_24_minus_23_is_calculated_before_sample_aggregation():
    recovery, joint = complete_population_rows()
    result = stage04.paired_comparisons(recovery, joint, [], CONFIG)
    row = next(
        x for x in result
        if x["focal_strand"] == "antisense"
        and x["metric"] == "varroa_2nt_reference_fraction_recovered"
    )
    # S1 median delta=(0.1,0.3)->0.2; S2 delta=0.8; across-sample median=0.5.
    assert math.isclose(row["sample_balanced_paired_delta_24_minus_23"], 0.5)


def test_spearman_known_vectors_and_ties():
    assert stage04.spearman_rho(list(range(16)), list(range(16))) == 1
    assert stage04.spearman_rho(list(range(16)), list(reversed(range(16)))) == -1
    assert stage04.spearman_rho([1, 1, 2], [1, 1, 2]) == 1


def test_redundancy_uses_matching_sixteen_features():
    sequence, general = [], []
    index = 0
    for position in stage04.POSITIONS:
        for nucleotide in stage04.NUCLEOTIDES:
            for mode in stage04.WEIGHTING_MODES:
                sequence.extend([
                    {
                        "focal_length": 23, "focal_strand": "antisense",
                        "weighting_mode": mode, "terminal_position": position,
                        "nucleotide": nucleotide, "metric": "E_joint_absolute",
                        "sample_balanced_median": index,
                    },
                    {
                        "focal_length": 23, "focal_strand": "antisense",
                        "weighting_mode": mode, "terminal_position": position,
                        "nucleotide": nucleotide,
                        "metric": "joint_vs_recovered_log2_contrast",
                        "sample_balanced_median": index,
                    },
                ])
                general.append({
                    "length": "23", "strand_scope": "antisense",
                    "weighting_mode": mode, "terminal_position": position,
                    "nucleotide": nucleotide,
                    "sample_balanced_median_enrichment_ratio": str(index),
                })
            index += 1
    result = stage04.calculate_redundancy(sequence, general)
    row = next(
        x for x in result
        if x["focal_length"] == 23 and x["focal_strand"] == "antisense"
        and x["comparison"] == "rho_joint_vs_general"
        and x["weighting_comparison"] == "abundance"
    )
    concordance = next(
        x for x in result
        if x["focal_length"] == 23 and x["focal_strand"] == "antisense"
        and x["comparison"] == "rho_joint_contrast_abundance_vs_unique"
    )
    assert row["n_matched_features"] == 16 and row["spearman_rho"] == 1
    assert concordance["n_matched_features"] == 16 and concordance["spearman_rho"] == 1
