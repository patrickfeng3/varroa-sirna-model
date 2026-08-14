import importlib.util
import math
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow/scripts/stage01.py"
SPEC = importlib.util.spec_from_file_location("stage01", SCRIPT)
stage01 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = stage01
SPEC.loader.exec_module(stage01)


CONFIG = stage01.Stage01Config(15, 35, 200, 20260814, "percentile", 0.95)


def eligibility(sample="S1", unit="V1", primary="true"):
    return {
        "sample": sample,
        "analysis_unit": unit,
        "biological_virus": unit,
        "polarity": "+ssRNA",
        "primary_eligible": primary,
    }


def feature(
    sample="S1", unit="V1", length=23, strand="sense", sequence="A" * 23,
    count=1, mapping_mode="exact", assignment="assigned",
):
    return {
        "sample": sample,
        "mapping_mode": mapping_mode,
        "virus": unit,
        "virus_assignment": assignment,
        "strand": strand,
        "sequence": sequence,
        "length": str(length),
        "count": str(count),
    }


def analyse(eligibility_rows=None, rows=None):
    return stage01.analyse_stage01(
        eligibility_rows or [eligibility()], rows or [], CONFIG
    )


def pair_length(result, mode, length, sample="S1", unit="V1"):
    return next(
        row for row in result["length_pair"]
        if row["sample"] == sample and row["analysis_unit"] == unit
        and row["weighting_mode"] == mode and row["length"] == length
    )


def pair_fraction(result, mode, sample="S1", unit="V1"):
    return next(
        row for row in result["fractions_pair"]
        if row["sample"] == sample and row["analysis_unit"] == unit
        and row["weighting_mode"] == mode
    )


def test_abundance_uses_count_not_row_count_and_unique_deduplicates():
    rows = [feature(count=5), feature(count=7)]
    result = analyse(rows=rows)
    assert pair_length(result, "abundance", 23)["length_count"] == 12
    assert pair_length(result, "unique_sequence", 23)["length_count"] == 1


def test_primary_eligibility_and_pair_matching_filter_rows():
    elig = [eligibility(), eligibility("S1", "V2", "false")]
    rows = [feature(count=3), feature(unit="V2", count=50), feature(unit="OTHER", count=100)]
    result = analyse(elig, rows)
    assert pair_length(result, "abundance", 23)["length_count"] == 3
    assert {row["analysis_unit"] for row in result["length_pair"]} == {"V1"}


def test_exact_assigned_and_valid_strand_filtering():
    rows = [
        feature(count=2),
        feature(count=20, mapping_mode="1mm"),
        feature(count=30, assignment="ambiguous_multi_virus"),
        feature(count=40, strand="ambiguous"),
    ]
    result = analyse(rows=rows)
    assert pair_length(result, "abundance", 23)["length_count"] == 2


def test_length_boundaries_are_inclusive_and_outside_rows_are_accounted():
    rows = [feature(length=14, sequence="A" * 14), feature(length=15, sequence="A" * 15),
            feature(length=35, sequence="C" * 35), feature(length=36, sequence="G" * 36)]
    result = analyse(rows=rows)
    assert pair_length(result, "abundance", 15)["length_count"] == 1
    assert pair_length(result, "abundance", 35)["length_count"] == 1
    outside = next(row for row in result["qc"] if row["metric"] == "rows_outside_15_35_nt")
    assert outside["value"] == 2


def test_unique_sequence_identity_is_strand_specific():
    seq = "A" * 23
    rows = [feature(sequence=seq, strand="sense"), feature(sequence=seq, strand="sense"),
            feature(sequence=seq, strand="antisense")]
    result = analyse(rows=rows)
    counts = next(row for row in result["counts_pair"] if row["weighting_mode"] == "unique_sequence")
    assert counts["n23_sense"] == 1
    assert counts["n23_antisense"] == 1
    assert counts["n23_total"] == 2


def test_identical_sequence_counts_in_different_analysis_units():
    elig = [eligibility(unit="V1"), eligibility(unit="V2")]
    rows = [feature(unit="V1"), feature(unit="V2")]
    result = analyse(elig, rows)
    assert pair_length(result, "unique_sequence", 23, unit="V1")["length_count"] == 1
    assert pair_length(result, "unique_sequence", 23, unit="V2")["length_count"] == 1


def test_length_fraction_uses_complete_15_35_denominator():
    result = analyse(rows=[feature(length=23, count=3), feature(length=24, count=1, sequence="C" * 24)])
    assert pair_length(result, "abundance", 23)["length_fraction"] == 0.75
    assert pair_length(result, "abundance", 24)["length_fraction"] == 0.25


def test_zero_denominator_produces_na_fraction_rank_and_indicators():
    result = analyse(rows=[])
    row = pair_length(result, "abundance", 23)
    assert row["length_fraction"] is None
    assert row["length_rank"] is None
    assert row["top1_indicator"] is None
    assert row["top3_indicator"] is None


def test_competition_ranking_ties_and_top_indicators():
    rows = [feature(length=23, count=5), feature(length=24, count=4, sequence="C" * 24),
            feature(length=25, count=4, sequence="G" * 25), feature(length=26, count=1, sequence="T" * 26)]
    result = analyse(rows=rows)
    assert pair_length(result, "abundance", 23)["length_rank"] == 1
    assert pair_length(result, "abundance", 24)["length_rank"] == 2
    assert pair_length(result, "abundance", 25)["length_rank"] == 2
    assert pair_length(result, "abundance", 26)["length_rank"] == 4
    assert pair_length(result, "abundance", 23)["top1_indicator"] == 1
    assert pair_length(result, "abundance", 24)["top3_indicator"] == 1
    assert pair_length(result, "abundance", 26)["top3_indicator"] == 0


def test_23_24_strand_fractions_delta_and_composition():
    rows = [
        feature(length=23, strand="sense", count=3, sequence="A" * 23),
        feature(length=23, strand="antisense", count=1, sequence="C" * 23),
        feature(length=24, strand="sense", count=2, sequence="G" * 24),
        feature(length=24, strand="antisense", count=6, sequence="T" * 24),
    ]
    row = pair_fraction(analyse(rows=rows), "abundance")
    assert row["sense_fraction_23"] == 0.75
    assert row["antisense_fraction_23"] == 0.25
    assert row["sense_fraction_24"] == 0.25
    assert row["antisense_fraction_24"] == 0.75
    assert row["delta_antisense_fraction_24_minus_23"] == 0.5
    assert row["length23_fraction_among_23_24"] == 1 / 3
    assert row["length24_fraction_among_23_24"] == 2 / 3


def test_23_24_zero_denominators_are_na_without_pseudocounts():
    row = pair_fraction(analyse(rows=[]), "abundance")
    for metric in stage01.FIXED_METRICS:
        assert row[metric] is None


def test_within_sample_median_across_viruses():
    elig = [eligibility(unit="V1"), eligibility(unit="V2")]
    rows = [
        feature(unit="V1", strand="sense", count=1), feature(unit="V1", strand="antisense", count=1, sequence="C" * 23),
        feature(unit="V2", strand="antisense", count=4, sequence="G" * 23),
    ]
    result = analyse(elig, rows)
    row = next(r for r in result["fixed_sample"] if r["weighting_mode"] == "abundance" and r["metric"] == "antisense_fraction_23")
    assert row["median_value"] == 0.75
    assert row["n_virus_units"] == 2


def test_sample_balanced_median_uses_sample_medians():
    elig = [eligibility("S1", "V1"), eligibility("S1", "V2"), eligibility("S2", "V1")]
    rows = [
        feature("S1", "V1", strand="sense", count=1),
        feature("S1", "V2", strand="antisense", count=1),
        feature("S2", "V1", strand="antisense", count=9), feature("S2", "V1", strand="sense", count=1, sequence="C" * 23),
    ]
    result = analyse(elig, rows)
    row = next(r for r in result["fixed_across"] if r["weighting_mode"] == "abundance" and r["metric"] == "antisense_fraction_23")
    assert row["sample_balanced_median"] == 0.7
    assert row["n_samples"] == 2
    assert row["n_sample_virus_units"] == 3


def test_sample_clustered_bootstrap_is_reproducible_with_fixed_seed():
    values = {"S1": 0.1, "S2": 0.5, "S3": 0.9}
    first = stage01.clustered_bootstrap_ci(values, 1000, 20260814, 0.95)
    second = stage01.clustered_bootstrap_ci(values, 1000, 20260814, 0.95)
    assert first == second
    assert first[2] == 1000
