from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/stage07.py"
SPEC = importlib.util.spec_from_file_location("stage07", SCRIPT)
assert SPEC and SPEC.loader
stage07 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage07)


def test_observed_antisense_is_not_reverse_complemented() -> None:
    sequence = "AACCGT"
    assert stage07.observed_physical_sequence(sequence, "antisense") == sequence
    assert stage07.observed_physical_sequence(sequence, "antisense") != stage07.reverse_complement(sequence)


def test_expected_antisense_is_reverse_complemented_and_sense_is_direct() -> None:
    window = "AACCGT"
    assert stage07.expected_physical_sequence(window, "sense") == window
    assert stage07.expected_physical_sequence(window, "antisense") == "ACGGTT"


def test_physical_position_and_three_prime_mapping() -> None:
    sequence = "ACGT"
    assert stage07.physical_position(sequence, 1) == "A"
    assert stage07.physical_position(sequence, 4) == "T"
    assert [stage07.position_from_3p(4, position) for position in range(1, 5)] == [4, 3, 2, 1]


def test_background_windows_never_cross_fasta_records() -> None:
    stats = stage07.enumerate_background(["A" * 13, "C" * 13], 14)
    assert stats["candidate"] == 0
    assert stats["valid"] == 0


def test_background_ambiguous_windows_are_excluded() -> None:
    stats = stage07.enumerate_background(["A" * 14, "A" * 13 + "N"], 14)
    assert stats["candidate"] == 2
    assert stats["valid"] == 1
    assert stats["excluded"] == 1


def test_repeated_background_windows_remain_separate_opportunities() -> None:
    stats = stage07.enumerate_background(["A" * 15], 14)
    assert stats["valid"] == 2
    assert stats["counts"]["sense"][(1, "A")] == 2
    assert stats["counts"]["antisense"][(1, "T")] == 2


def test_representation_fractions_enrichment_and_delta() -> None:
    observed = stage07.safe_fraction(3, 4)
    expected = stage07.safe_fraction(1, 2)
    assert observed == 0.75
    assert expected == 0.5
    assert stage07.representation_enrichment(observed, expected, 4) == 1.5
    assert stage07.safe_delta(observed, expected) == 0.25


def test_representation_zero_denominator_is_na_without_pseudocount() -> None:
    assert stage07.safe_fraction(1, 0) is None
    assert stage07.representation_enrichment(0.5, 0.0, 4) is None
    assert stage07.representation_enrichment(None, 0.5, 0) is None


def test_accumulation_ratio_delta_and_zero_denominator() -> None:
    ratio, delta, log_ratio = stage07.accumulation_metrics(0.25, 0.5)
    assert ratio == 2.0
    assert delta == 0.25
    assert log_ratio == 1.0
    ratio, delta, log_ratio = stage07.accumulation_metrics(0.0, 0.5)
    assert ratio is None
    assert delta == 0.5
    assert log_ratio is None


def test_continuous_gc9_14() -> None:
    sequence = "A" * 8 + "GCGCAA" + "T" * 10
    assert stage07.gc9_14_fraction(sequence) == pytest.approx(4 / 6)


def test_exact_a10_extraction_and_literature_row() -> None:
    sequence = "C" * 9 + "A" + "G" * 13
    assert stage07.physical_position(sequence, 10) == "A"
    positional = [{
        "strand": "antisense", "length": 23, "position_5p": 10, "nucleotide": "A",
        "weighting_mode": "unique_sequence", "endpoint": "unique_representation",
        "sample_balanced_representation_enrichment": 1.2,
        "sample_balanced_representation_delta_fraction": 0.1,
        "bootstrap_ci_low": 1.0, "bootstrap_ci_high": 1.4, "n_samples_total": 20,
        "sign_test_n_nonzero": 18, "sign_test_n_positive": 12,
        "sign_test_n_negative": 6, "sign_test_estimability": "estimable", "raw_p": 0.2,
    }]
    literature = stage07.build_literature_validation(positional, [])
    assert len(literature) == 1
    assert literature[0]["feature"] == "A10"
    assert literature[0]["effect_estimate"] == 1.2


def pair_row(
    sample: str,
    enrichment: float,
    delta: float,
    analysis_unit: str,
) -> dict[str, object]:
    return {
        "sample": sample, "analysis_unit": analysis_unit, "biological_virus": analysis_unit,
        "polarity": "positive", "length": 23, "strand": "antisense",
        "weighting_mode": "unique_sequence", "position_5p": 10,
        "position_from_3p": 14, "nucleotide": "A", "observed_nucleotide_weight": 1.0,
        "observed_total_weight": 2.0, "observed_fraction": 0.5,
        "expected_nucleotide_weight": 1.0, "expected_total_windows": 4,
        "expected_fraction": 0.25, "representation_enrichment": enrichment,
        "representation_delta_fraction": delta, "unique_fraction": 0.5,
        "abundance_fraction": 0.6, "accumulation_ratio": 1.2,
        "accumulation_delta_fraction": 0.1, "log2_accumulation_ratio": math.log2(1.2),
        "valid_background_windows": 4,
    }


def test_pair_to_sample_to_dataset_median_hierarchy() -> None:
    rows = [
        pair_row("sample_a", 2.0, 0.1, "v1"),
        pair_row("sample_a", 4.0, 0.3, "v2"),
        pair_row("sample_b", 10.0, 0.5, "v1"),
    ]
    samples, summary = stage07.aggregate_positional(rows, bootstrap_replicates=20, seed=7)
    assert len(samples) == 2
    assert summary[0]["sample_balanced_representation_enrichment"] == 6.5
    assert summary[0]["sample_balanced_representation_delta_fraction"] == pytest.approx(0.35)


def test_sample_clustered_bootstrap_is_reproducible() -> None:
    values = {"a": 1.0, "b": 2.0, "c": 10.0}
    first = stage07.sample_clustered_bootstrap(values, replicates=100, seed=42)
    second = stage07.sample_clustered_bootstrap(values, replicates=100, seed=42)
    assert first == second
    assert first[2] == 100


def test_two_sided_exact_sign_test() -> None:
    result = stage07.exact_sign_test([1, 2, 3, -1])
    assert result["n_nonzero"] == 4
    assert result["n_positive"] == 3
    assert result["n_negative"] == 1
    assert result["raw_p"] == 0.625
    empty = stage07.exact_sign_test([0, 0])
    assert empty["raw_p"] is None
    assert empty["estimability"] == "not_estimable_no_nonzero_sample_deltas"


def test_bh_correction() -> None:
    adjusted = stage07.adjust_pvalues([0.01, 0.04, 0.03], "BH")
    assert adjusted == pytest.approx([0.03, 0.04, 0.04])


def test_by_correction_is_more_conservative_than_bh() -> None:
    raw = [0.01, 0.04, 0.03]
    bh = stage07.adjust_pvalues(raw, "BH")
    by = stage07.adjust_pvalues(raw, "BY")
    assert all(by_value >= bh_value for by_value, bh_value in zip(by, bh))


def test_discovery_corrections_are_separate_by_family() -> None:
    rows = [
        {"length": 23, "endpoint": "unique_representation", "strand": "antisense", "position_5p": 3, "raw_p": 0.01, "bh_p": None, "by_p": None},
        {"length": 23, "endpoint": "unique_representation", "strand": "antisense", "position_5p": 4, "raw_p": 0.04, "bh_p": None, "by_p": None},
        {"length": 24, "endpoint": "unique_representation", "strand": "antisense", "position_5p": 3, "raw_p": 0.03, "bh_p": None, "by_p": None},
        {"length": 23, "endpoint": "accumulation", "strand": "sense", "position_5p": 3, "raw_p": 0.02, "bh_p": None, "by_p": None},
    ]
    stage07.apply_discovery_corrections(rows)
    assert rows[0]["bh_p"] == 0.02
    assert rows[1]["bh_p"] == 0.04
    assert rows[2]["bh_p"] == 0.03
    assert rows[3]["bh_p"] == 0.02


def valid_scope() -> tuple[list[dict[str, str]], dict[tuple[str, str], list[str]]]:
    eligibility: list[dict[str, str]] = []
    pairs: list[tuple[str, str]] = []
    for sample_index in range(20):
        units = 3 if sample_index < 14 else 2
        for unit_index in range(units):
            pair = (f"sample_{sample_index:02d}", f"virus_{unit_index}")
            pairs.append(pair)
            eligibility.append({
                "sample": pair[0], "analysis_unit": pair[1], "primary_eligible": "true",
                "biological_virus": pair[1], "polarity": "positive",
            })
    assert len(pairs) == 54
    backgrounds = {pair: ["ACGT" * 8] for pair in pairs}
    return eligibility, backgrounds


def test_unique_collapse_and_abundance_count_accounting() -> None:
    eligibility, backgrounds = valid_scope()
    pair = (eligibility[0]["sample"], eligibility[0]["analysis_unit"])
    sequence = "A" + "C" * 22
    feature_rows = [
        {"sample": pair[0], "virus": pair[1], "mapping_mode": "exact", "virus_assignment": "assigned", "strand": "antisense", "length": "23", "sequence": sequence, "count": "2"},
        {"sample": pair[0], "virus": pair[1], "mapping_mode": "exact", "virus_assignment": "assigned", "strand": "antisense", "length": "23", "sequence": sequence, "count": "3"},
    ]
    rows, _gc, _stats = stage07.build_pair_tables(eligibility, feature_rows, backgrounds)
    selected = [
        row for row in rows
        if row["sample"] == pair[0] and row["analysis_unit"] == pair[1]
        and row["length"] == 23 and row["strand"] == "antisense"
        and row["position_5p"] == 1 and row["nucleotide"] == "A"
    ]
    by_mode = {row["weighting_mode"]: row for row in selected}
    assert by_mode["unique_sequence"]["observed_total_weight"] == 1.0
    assert by_mode["abundance"]["observed_total_weight"] == 5.0
    assert by_mode["unique_sequence"]["observed_fraction"] == 1.0
    assert by_mode["abundance"]["observed_fraction"] == 1.0


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_stage02_four_terminal_regression_within_tolerance(tmp_path: Path) -> None:
    terminal_positions = {"5p1": 1, "5p2": 2, "3p2": 22, "3p1": 23}
    pair_rows = []
    summary_rows = []
    old_pair = []
    old_across = []
    for terminal, position in terminal_positions.items():
        pair_rows.append({
            "sample": "s", "analysis_unit": "v", "length": 23, "strand": "antisense",
            "weighting_mode": "unique_sequence", "position_5p": position, "nucleotide": "A",
            "observed_fraction": 0.5, "expected_fraction": 0.25,
            "representation_enrichment": 2.0,
        })
        summary_rows.append({
            "length": 23, "strand": "antisense", "weighting_mode": "unique_sequence",
            "position_5p": position, "nucleotide": "A",
            "sample_balanced_representation_enrichment": 2.0,
        })
        old_pair.append({
            "sample": "s", "analysis_unit": "v", "length": 23,
            "strand_scope": "antisense", "weighting_mode": "unique_sequence",
            "terminal_position": terminal, "nucleotide": "A", "observed_fraction": 0.5,
            "expected_fraction": 0.25, "enrichment_ratio": 2.0,
        })
        old_across.append({
            "length": 23, "strand_scope": "antisense", "weighting_mode": "unique_sequence",
            "terminal_position": terminal, "nucleotide": "A",
            "sample_balanced_median_enrichment_ratio": 2.0,
        })
    pair_path = tmp_path / "pair.tsv"
    across_path = tmp_path / "across.tsv"
    write_tsv(pair_path, old_pair)
    write_tsv(across_path, old_across)
    before_pair = pair_path.read_bytes()
    before_across = across_path.read_bytes()
    identities = stage07.validate_stage02_references(pair_path, across_path)
    checks, passed, maximum, count = stage07.stage02_terminal_regression(
        pair_rows, summary_rows, pair_path, across_path
    )
    assert passed
    assert count == 16
    assert maximum <= 1e-12
    assert {row["position_5p"] for row in checks} == {1, 2, 22, 23}
    assert pair_path.read_bytes() == before_pair
    assert across_path.read_bytes() == before_across
    assert len(identities["stage02_pair_reference_sha256"]) == 64
    assert len(identities["stage02_across_reference_sha256"]) == 64


def test_missing_stage02_regression_reference_fails_clearly(tmp_path: Path) -> None:
    existing = tmp_path / "existing.tsv"
    existing.write_text("header\n", encoding="utf-8")
    missing = tmp_path / "missing.tsv"
    with pytest.raises(stage07.Stage07Error, match="required Stage 02 regression reference is missing"):
        stage07.validate_stage02_references(existing, missing)


def test_stage07_rule_keeps_stage02_references_out_of_dag_inputs() -> None:
    rule_text = (ROOT / "workflow/rules/stage07.smk").read_text(encoding="utf-8")
    input_block = rule_text.split("    output:", 1)[0]
    assert "terminal_enrichment_by_pair.tsv" not in input_block
    assert "terminal_enrichment_across_dataset.tsv" not in input_block
    assert "stage02_pair_reference=" in rule_text
    assert "stage02_across_reference=" in rule_text


@pytest.mark.parametrize(("length", "expected"), [(23, 18), (24, 19)])
def test_regional_gc6_window_counts(length: int, expected: int) -> None:
    windows = stage07.regional_windows(length)
    assert len(windows) == expected
    assert windows[0]["start_5p"] == 1
    assert windows[-1]["end_5p"] == length


def test_regional_gc6_coordinates_and_labels() -> None:
    window = stage07.regional_windows(23)[13]
    assert window == {
        "start_5p": 14,
        "end_5p": 19,
        "near_3p": 5,
        "far_3p": 10,
        "region_5p": "GC14-19",
        "region_3p": "GC_3p5-10",
    }


def test_regional_gc_rejects_noncanonical_widths() -> None:
    with pytest.raises(ValueError, match="width must be 6"):
        stage07.regional_windows(23, width=5)


def test_direct_regional_gc6_fraction() -> None:
    sequence = "A" * 8 + "GCGCAA" + "T" * 9
    assert stage07.regional_gc6_fraction(sequence, 9) == pytest.approx(4 / 6)


def test_positional_gc_derivation_equals_direct_sequence_mean() -> None:
    sequences = [
        "AACCGGTTAACCGGTTAACCGGT",
        "GGCCAATTGGCCAATTGGCCAAT",
        "TTGGCCAATTGGCCAATTGGCCA",
    ]
    start = 9
    direct = sum(stage07.regional_gc6_fraction(sequence, start) for sequence in sequences) / len(sequences)
    positional = {
        position: sum(sequence[position - 1] in {"G", "C"} for sequence in sequences) / len(sequences)
        for position in range(1, 24)
    }
    derived = stage07.regional_gc_from_positional_fractions(positional, start)
    assert derived == pytest.approx(direct, abs=1e-12)


def synthetic_positional_rows(length: int = 23) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    distributions = {
        "unique_sequence": {"A": 0.1, "C": 0.2, "G": 0.3, "T": 0.4},
        "abundance": {"A": 0.2, "C": 0.3, "G": 0.1, "T": 0.4},
    }
    for mode, distribution in distributions.items():
        for position in range(1, length + 1):
            for nucleotide, observed in distribution.items():
                rows.append({
                    "sample": "sample", "analysis_unit": "virus",
                    "biological_virus": "virus", "polarity": "positive",
                    "length": length, "strand": "antisense", "weighting_mode": mode,
                    "position_5p": position, "nucleotide": nucleotide,
                    "observed_fraction": observed, "expected_fraction": 0.25,
                })
    return rows


def test_regional_pair_deltas_derived_from_existing_positional_rows() -> None:
    regional = stage07.build_regional_gc_pair_rows(synthetic_positional_rows())
    assert len(regional) == 18
    first = regional[0]
    assert first["observed_gc6_mean_unique"] == pytest.approx(0.5)
    assert first["observed_gc6_mean_abundance"] == pytest.approx(0.4)
    assert first["expected_gc6_mean"] == pytest.approx(0.5)
    assert first["regional_gc6_delta_unique_vs_expected"] == pytest.approx(0.0)
    assert first["regional_gc6_delta_abundance_vs_expected"] == pytest.approx(-0.1)
    assert first["regional_gc6_accumulation_delta"] == pytest.approx(-0.1)


def regional_pair_fixture(
    sample: str,
    unit: str,
    start: int,
    unique: float,
    abundance: float,
    expected: float,
    length: int = 23,
    strand: str = "antisense",
) -> dict[str, object]:
    coordinates = stage07.regional_windows(length)[start - 1]
    unique_delta, abundance_delta, accumulation_delta = stage07.regional_gc_deltas(
        unique, abundance, expected
    )
    return {
        "sample": sample, "analysis_unit": unit, "biological_virus": unit,
        "polarity": "positive", "strand": strand, "length": length, **coordinates,
        "observed_gc6_mean_unique": unique,
        "observed_gc6_mean_abundance": abundance,
        "expected_gc6_mean": expected,
        "regional_gc6_delta_unique_vs_expected": unique_delta,
        "regional_gc6_delta_abundance_vs_expected": abundance_delta,
        "regional_gc6_accumulation_delta": accumulation_delta,
    }


def test_regional_pair_sample_dataset_uses_existing_median_hierarchy() -> None:
    rows = [
        regional_pair_fixture("sample_a", "v1", 1, 0.5, 0.5, 0.4),
        regional_pair_fixture("sample_a", "v2", 1, 0.7, 0.7, 0.4),
        regional_pair_fixture("sample_b", "v1", 1, 0.9, 0.9, 0.4),
    ]
    samples, summaries, _discovery = stage07.aggregate_regional_gc(
        rows, [], bootstrap_replicates=20, seed=7
    )
    unique = next(row for row in summaries if row["endpoint"] == "unique_representation")
    assert len(samples) == 2
    assert unique["sample_balanced_regional_gc6_delta"] == pytest.approx(0.35)


def test_regional_gc9_14_regresses_to_existing_gc_engine() -> None:
    regional_pair = [
        regional_pair_fixture("sample_a", "v1", 9, 0.6, 0.5, 0.4),
        regional_pair_fixture("sample_a", "v2", 9, 0.5, 0.45, 0.4),
        regional_pair_fixture("sample_b", "v1", 9, 0.7, 0.6, 0.4),
    ]
    gc_pair = [{
        "sample": row["sample"], "analysis_unit": row["analysis_unit"],
        "biological_virus": row["biological_virus"], "polarity": row["polarity"],
        "length": row["length"], "strand": row["strand"],
        "observed_GC9_14_mean_unique": row["observed_gc6_mean_unique"],
        "observed_GC9_14_mean_abundance": row["observed_gc6_mean_abundance"],
        "expected_GC9_14_mean": row["expected_gc6_mean"],
        "GC9_14_delta_unique_vs_expected": row["regional_gc6_delta_unique_vs_expected"],
        "GC9_14_delta_abundance_vs_expected": row["regional_gc6_delta_abundance_vs_expected"],
        "GC9_14_accumulation_delta": row["regional_gc6_accumulation_delta"],
    } for row in regional_pair]
    gc_sample, gc_summary = stage07.aggregate_gc(gc_pair, bootstrap_replicates=20, seed=9)
    regional_sample, regional_summary, _ = stage07.aggregate_regional_gc(
        regional_pair, [], bootstrap_replicates=20, seed=9
    )
    passed, maximum, count = stage07.regional_gc9_14_regression(
        regional_pair, regional_sample, regional_summary, gc_pair, gc_sample, gc_summary
    )
    assert passed
    assert maximum <= 1e-12
    assert count == 33


def test_gc9_14_excluded_and_regional_family_sizes_are_exact() -> None:
    rows: list[dict[str, object]] = []
    for length in (23, 24):
        for endpoint in ("unique_representation", "abundance_representation", "accumulation"):
            for strand in ("antisense", "sense"):
                for coordinates in stage07.regional_windows(length):
                    start = coordinates["start_5p"]
                    rows.append({
                        "length": length, "endpoint": endpoint, "strand": strand,
                        **coordinates, "raw_p": 0.01 + int(start) / 1000,
                        "evidence_class": (
                            "literature_validation_gc9_14"
                            if start == 9 else "exploratory_regional_gc6"
                        ),
                        "regional_bh_p": None, "regional_by_p": None,
                    })
    stage07.apply_regional_gc_corrections(rows)
    gc9 = [row for row in rows if row["start_5p"] == 9]
    assert all(row["regional_bh_p"] is None and row["regional_by_p"] is None for row in gc9)
    exploratory = [row for row in rows if row["evidence_class"] == "exploratory_regional_gc6"]
    for length, expected in ((23, 17), (24, 18)):
        families = {}
        for row in exploratory:
            if row["length"] == length:
                families.setdefault((row["endpoint"], row["strand"]), []).append(row)
        assert len(families) == 6
        assert {len(family) for family in families.values()} == {expected}
    assert all(row["regional_bh_p"] is not None and row["regional_by_p"] is not None for row in exploratory)
