import importlib.util
import math
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow/scripts/stage02.py"
SPEC = importlib.util.spec_from_file_location("stage02", SCRIPT)
stage02 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = stage02
SPEC.loader.exec_module(stage02)


CONFIG = stage02.Stage02Config(
    (23, 24), ("5p1", "5p2", "3p2", "3p1"), 10, 20260814,
    "percentile", 0.95, 1e-12,
)


def eligibility(sample="S1", unit="V1"):
    return {
        "sample": sample, "analysis_unit": unit, "biological_virus": unit,
        "polarity": "+ssRNA", "primary_eligible": "true",
    }


def feature(
    sample="S1", unit="V1", length=23, strand="sense", sequence=None,
    count=1, mapping_mode="exact", assignment="assigned",
):
    sequence = sequence if sequence is not None else "A" * length
    return {
        "sample": sample, "virus": unit, "mapping_mode": mapping_mode,
        "virus_assignment": assignment, "strand": strand,
        "length": str(length), "sequence": sequence, "count": str(count),
    }


def analyse(rows, eligibility_rows=None, backgrounds=None):
    eligibility_rows = eligibility_rows or [eligibility()]
    backgrounds = backgrounds or {("S1", "V1"): ["A" * 24]}
    return stage02.analyse_stage02(eligibility_rows, rows, backgrounds, CONFIG)


def find_row(rows, **wanted):
    return next(row for row in rows if all(row[key] == value for key, value in wanted.items()))


def test_terminal_coordinate_extraction():
    assert stage02.terminal_bases("ACAAAAATG") == {"5p1": "A", "5p2": "C", "3p2": "T", "3p1": "G"}


def test_observed_antisense_is_not_reverse_complemented():
    sequence = "AC" + "A" * 19 + "TG"
    result = analyse([feature(strand="antisense", sequence=sequence)])
    for position, nucleotide in {"5p1": "A", "5p2": "C", "3p2": "T", "3p1": "G"}.items():
        row = find_row(result["observed"], strand_scope="antisense", weighting_mode="abundance", length=23, terminal_position=position, nucleotide=nucleotide)
        assert row["observed_fraction"] == 1


def test_expected_antisense_is_reverse_complemented():
    stats = stage02.enumerate_background(["AAAC"], 4)
    assert stats["valid"] == 1
    expected = stage02.terminal_bases("GTTT")
    for position, nucleotide in expected.items():
        assert stats["counts"]["antisense"][(position, nucleotide)] == 1


def test_abundance_uses_count_not_row_count():
    result = analyse([feature(count=7), feature(count=13, sequence="C" * 23)])
    row = find_row(result["observed"], strand_scope="sense", weighting_mode="abundance", length=23, terminal_position="5p1", nucleotide="A")
    assert row["observed_total_weight"] == 20
    assert row["observed_terminal_weight"] == 7


def test_unique_sequence_deduplication_is_strand_specific():
    sequence = "A" * 23
    result = analyse([feature(sequence=sequence), feature(sequence=sequence), feature(sequence=sequence, strand="antisense")])
    sense = find_row(result["observed"], strand_scope="sense", weighting_mode="unique_sequence", length=23, terminal_position="5p1", nucleotide="A")
    antisense = find_row(result["observed"], strand_scope="antisense", weighting_mode="unique_sequence", length=23, terminal_position="5p1", nucleotide="A")
    combined = find_row(result["observed"], strand_scope="combined", weighting_mode="unique_sequence", length=23, terminal_position="5p1", nucleotide="A")
    assert sense["observed_total_weight"] == 1
    assert antisense["observed_total_weight"] == 1
    assert combined["observed_total_weight"] == 2


def test_repeated_identical_background_windows_are_separate_opportunities():
    stats = stage02.enumerate_background(["AAAA"], 3)
    assert stats["valid"] == 2
    assert stats["counts"]["sense"][("5p1", "A")] == 2


def test_background_windows_never_cross_records():
    stats = stage02.enumerate_background(["AA", "AA"], 3)
    assert stats["candidate"] == 0
    assert stats["valid"] == 0


def test_background_windows_containing_n_are_excluded():
    stats = stage02.enumerate_background(["AANA"], 3)
    assert stats["candidate"] == 2
    assert stats["valid"] == 0
    assert stats["excluded"] == 2


def test_observed_and_expected_fractions_sum_to_one():
    result = analyse([feature(sequence="A" * 23), feature(sequence="C" * 23, count=2)])
    for table, field in ((result["observed"], "observed_fraction"), (result["expected"], "expected_fraction")):
        groups = {}
        for row in table:
            if row[field] is None:
                continue
            key = tuple(row[name] for name in ("sample", "analysis_unit", "length", "strand_scope", "weighting_mode", "terminal_position"))
            groups.setdefault(key, []).append(row[field])
        assert all(math.isclose(sum(values), 1.0, abs_tol=1e-12) for values in groups.values())


def test_combined_expectation_uses_observed_strand_weights_not_half():
    rows = [feature(strand="sense", count=3), feature(strand="antisense", count=1, sequence="C" * 23)]
    result = analyse(rows)
    expected_a = find_row(result["expected"], strand_scope="combined", weighting_mode="abundance", length=23, terminal_position="5p1", nucleotide="A")
    expected_t = find_row(result["expected"], strand_scope="combined", weighting_mode="abundance", length=23, terminal_position="5p1", nucleotide="T")
    assert expected_a["expected_fraction"] == 0.75
    assert expected_t["expected_fraction"] == 0.25


def test_abundance_and_unique_modes_can_have_different_combined_expectations():
    rows = [
        feature(strand="sense", count=100, sequence="A" * 23),
        feature(strand="antisense", count=1, sequence="C" * 23),
        feature(strand="antisense", count=1, sequence="G" * 23),
    ]
    result = analyse(rows)
    abundance = find_row(result["expected"], strand_scope="combined", weighting_mode="abundance", length=23, terminal_position="5p1", nucleotide="A")
    unique = find_row(result["expected"], strand_scope="combined", weighting_mode="unique_sequence", length=23, terminal_position="5p1", nucleotide="A")
    assert math.isclose(abundance["expected_fraction"], 100 / 102)
    assert math.isclose(unique["expected_fraction"], 1 / 3)


def test_enrichment_edge_cases_without_pseudocounts():
    assert stage02.enrichment_ratio(0.0, 0.25, 10) == 0
    assert stage02.enrichment_ratio(0.5, 0.0, 10) is None
    assert stage02.enrichment_ratio(None, 0.25, 0) is None


def test_sample_balanced_median_is_virus_then_sample_median():
    rows = []
    for sample, unit, value in (("S1", "V1", 0.0), ("S1", "V2", 1.0), ("S2", "V1", 0.9)):
        rows.append({
            "sample": sample, "analysis_unit": unit, "length": 23,
            "strand_scope": "antisense", "weighting_mode": "abundance",
            "terminal_position": "5p1", "nucleotide": "A", "enrichment_ratio": value,
        })
    sample_rows, across = stage02.aggregate_enrichment(rows, CONFIG)
    s1 = find_row(sample_rows, sample="S1")
    assert s1["sample_enrichment_median"] == 0.5
    assert across[0]["sample_balanced_median_enrichment_ratio"] == 0.7


def test_sample_clustered_bootstrap_is_reproducible():
    values = {"S1": [0.1, 0.9], "S2": [0.4], "S3": [0.8]}
    first = stage02.sample_clustered_bootstrap(values, 1000, 20260814, 0.95)
    second = stage02.sample_clustered_bootstrap(values, 1000, 20260814, 0.95)
    assert first == second
    assert first[2] == 1000


def test_pooled_expected_fraction_is_abundance_weighted():
    rows = [
        {"length": 23, "strand_scope": "antisense", "weighting_mode": "abundance", "terminal_position": "5p1", "nucleotide": "A", "sample": "S1", "analysis_unit": "V1", "observed_total_weight": 10.0, "observed_terminal_weight": 1.0, "expected_fraction": 0.2},
        {"length": 23, "strand_scope": "antisense", "weighting_mode": "abundance", "terminal_position": "5p1", "nucleotide": "A", "sample": "S2", "analysis_unit": "V1", "observed_total_weight": 90.0, "observed_terminal_weight": 45.0, "expected_fraction": 0.8},
    ]
    result = stage02.pooled_abundance(rows)[0]
    assert result["pooled_abundance_observed_fraction"] == 0.46
    assert result["pooled_abundance_expected_fraction"] == 0.74
    assert math.isclose(result["pooled_abundance_enrichment_ratio"], 0.46 / 0.74)


def test_known_synthetic_spearman_vectors():
    ascending = list(range(16))
    assert stage02.spearman_rho(ascending, ascending) == 1
    assert stage02.spearman_rho(ascending, list(reversed(ascending))) == -1
