"""Targeted deterministic tests for canonical Stage 08 only."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/stage08.py"
RESOURCE = ROOT / "resources/parameters/zuber_2022_wcf_dg37.tsv"
RULE = ROOT / "workflow/rules/stage08.smk"
SPEC = importlib.util.spec_from_file_location("stage08", SCRIPT)
assert SPEC and SPEC.loader
STAGE08 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(STAGE08)


EXPECTED_STACKS = {
    "AA": -0.94, "AC": -2.25, "AG": -2.01, "AU": -1.09,
    "CA": -2.07, "CC": -3.28, "CG": -2.33, "CU": -2.01,
    "GA": -2.42, "GC": -3.46, "GG": -3.28, "GU": -2.25,
    "UA": -1.29, "UC": -2.42, "UG": -2.07, "UU": -0.94,
}
EXPECTED_CORRECTIONS = {
    "AU_terminal_on_AU_penultimate": 0.22,
    "AU_terminal_on_GC_penultimate": 0.44,
}


def _write_tsv(path: Path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _synthetic_registry(tmp_path: Path):
    sequence = "ACGTACGTAC"
    fasta = tmp_path / "transcript.fa"
    fasta.write_text(">record_a\n" + sequence + "\n", encoding="utf-8")
    manifest = tmp_path / "targets.tsv"
    fields = [
        "target_id", "transcript_id", "display_name", "organism", "molecule_type",
        "fasta_path", "fasta_record_id", "annotation_path", "expected_length_nt",
        "sequence_sha256_uppercase_dna", "candidate_lengths_nt", "source_database",
        "source_accession_version",
    ]
    _write_tsv(manifest, fields, [{
        "target_id": "generic_target", "transcript_id": "tx_any",
        "display_name": "Synthetic", "organism": "Example", "molecule_type": "mRNA",
        "fasta_path": str(fasta), "fasta_record_id": "record_a", "annotation_path": "NA",
        "expected_length_nt": len(sequence),
        "sequence_sha256_uppercase_dna": hashlib.sha256(sequence.encode()).hexdigest(),
        "candidate_lengths_nt": "5", "source_database": "synthetic",
        "source_accession_version": "tx_any.1",
    }])
    candidates = tmp_path / "candidates.tsv"
    rows = []
    for index, start in enumerate((2, 3), start=1):
        end = start + 4
        target_dna = sequence[start - 1:end]
        target_rna = target_dna.replace("T", "U")
        rows.append({
            "target_id": "generic_target", "transcript_id": "tx_any",
            "display_name": "Synthetic", "organism": "Example",
            "candidate_id": f"generic_{index}", "candidate_length_nt": 5,
            "start_1based": start, "end_1based": end,
            "target_sequence_dna": target_dna, "target_sequence_rna": target_rna,
            "antisense_guide_sequence_rna": STAGE08.reverse_complement_rna(target_rna),
            "annotation_status": "unavailable", "start_region": "NA", "end_region": "NA",
            "overlap_regions": "NA", "crosses_annotation_boundary": "NA",
        })
    _write_tsv(candidates, STAGE08.STAGE06_COLUMNS, rows)
    return manifest, candidates, rows


def test_generic_candidate_slice_orientation_and_count_are_preserved(tmp_path):
    manifest, candidate_path, source_rows = _synthetic_registry(tmp_path)
    transcripts = STAGE08.load_registered_transcripts(manifest)
    candidates = STAGE08.load_and_validate_candidates(candidate_path, transcripts)
    assert len(candidates) == len(source_rows) == 2
    assert candidates[0]["target_sequence_dna"] == "CGTAC"
    assert candidates[0]["antisense_guide_sequence_rna"] == "GUACG"
    assert candidates[0]["target_id"] == "generic_target"


def test_candidate_target_slice_mismatch_is_rejected(tmp_path):
    manifest, candidate_path, rows = _synthetic_registry(tmp_path)
    rows[0]["target_sequence_dna"] = "AAAAA"
    rows[0]["target_sequence_rna"] = "AAAAA"
    rows[0]["antisense_guide_sequence_rna"] = "UUUUU"
    _write_tsv(candidate_path, STAGE08.STAGE06_COLUMNS, rows)
    with pytest.raises(STAGE08.Stage08Error, match="transcript slice mismatch"):
        STAGE08.load_and_validate_candidates(
            candidate_path, STAGE08.load_registered_transcripts(manifest)
        )


def test_candidate_guide_reverse_complement_mismatch_is_rejected(tmp_path):
    manifest, candidate_path, rows = _synthetic_registry(tmp_path)
    rows[0]["antisense_guide_sequence_rna"] = "AAAAA"
    _write_tsv(candidate_path, STAGE08.STAGE06_COLUMNS, rows)
    with pytest.raises(STAGE08.Stage08Error, match="reverse-complement mismatch"):
        STAGE08.load_and_validate_candidates(
            candidate_path, STAGE08.load_registered_transcripts(manifest)
        )


def test_lunp_whole_and_seed_interval_indexing(tmp_path):
    lunp = tmp_path / "synthetic_lunp"
    lines = ["# i followed by probabilities for lengths 1..24"]
    for end in range(1, 31):
        values = [f"{end / 1000 + length / 10000:.6f}" for length in range(1, 25)]
        lines.append(str(end) + "\t" + "\t".join(values))
    lunp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    table = STAGE08.parse_lunp(lunp)
    assert table.interval_probability(2, 5) == pytest.approx(0.0054)
    # Candidate end e=15 gives target g2-g8 interval [8,14]: row/end 14, length 7.
    assert table.interval_probability(8, 14) == pytest.approx(0.0147)
    assert table.interval_probability(8, 14) != table.interval_probability(9, 15)


@pytest.mark.parametrize("value", [-0.01, 1.01, float("inf"), float("nan")])
def test_accessibility_probability_validation(value):
    with pytest.raises(STAGE08.Stage08Error, match="invalid accessibility"):
        STAGE08.validate_probability(value)


def test_zero_accessibility_is_valid_without_pseudocount():
    table = STAGE08.LunpTable({(8, 7): 0.0})
    assert table.interval_probability(2, 8) == 0.0


@pytest.mark.parametrize(
    ("requested_w", "requested_l", "length", "expected"),
    [(150, 100, 710, (150, 100)), (150, 100, 90, (90, 89)), (100, 80, 50, (50, 49))],
)
def test_effective_window_rule(requested_w, requested_l, length, expected):
    assert STAGE08.effective_window_parameters(requested_w, requested_l, length) == expected


def test_one_transcript_fold_is_reused_across_candidates_and_seed_uses_same_table():
    candidates = [
        {"target_id": "t", "transcript_id": "tx", "candidate_id": "c1",
         "candidate_length_nt": 8, "start_1based": 1, "end_1based": 8},
        {"target_id": "t", "transcript_id": "tx", "candidate_id": "c2",
         "candidate_length_nt": 8, "start_1based": 2, "end_1based": 9},
    ]
    transcripts = {("t", "tx"): {"sequence_rna": "ACGUACGUACGU"}}
    calls = []

    def fake_runner(sequence, window, span, ulength, temperature):
        calls.append((sequence, window, span, ulength, temperature))
        return STAGE08.LunpTable({(8, 8): 0.1, (9, 8): 0.2, (7, 7): 0.3, (8, 7): 0.4})

    parameters = [
        {"id": name, "window_nt": window, "max_bp_span_nt": span}
        for name, window, span in (
            ("W150_L100_main", 150, 100),
            ("W100_L80_sensitivity", 100, 80),
            ("W200_L150_sensitivity", 200, 150),
        )
    ]
    values, _, run_count = STAGE08.calculate_accessibilities(
        candidates, transcripts, parameters, 37.0, fake_runner
    )
    assert run_count == len(calls) == 3
    assert values["c1"]["target_whole_p_unpaired"] == 0.1
    assert values["c1"]["target_seed_g2_8_p_unpaired"] == 0.3
    assert values["c2"]["target_seed_g2_8_p_unpaired"] == 0.4


def test_zuber_resource_contains_exact_constants_and_provenance():
    stacks, corrections = STAGE08.load_zuber_parameters(RESOURCE)
    assert stacks == EXPECTED_STACKS
    assert corrections == EXPECTED_CORRECTIONS
    text = RESOURCE.read_text(encoding="utf-8")
    assert "10.1093/nar/gkac261" in text
    assert "kcal/mol" in text


def test_reverse_complement_equivalent_zuber_stacks_agree():
    stacks, _ = STAGE08.load_zuber_parameters(RESOURCE)
    for dinucleotide, value in stacks.items():
        assert stacks[STAGE08.reverse_complement_rna(dinucleotide)] == value


def test_four_and_five_bp_use_exact_stack_counts_without_initiation_or_symmetry():
    stacks = {left + right: 1.0 for left in "ACGU" for right in "ACGU"}
    corrections = {
        "AU_terminal_on_AU_penultimate": 0.0,
        "AU_terminal_on_GC_penultimate": 0.0,
    }
    assert STAGE08.terminal_dg("AAAAA", 4, stacks, corrections) == 3.0
    assert STAGE08.terminal_dg("AAAAA", 5, stacks, corrections) == 4.0


def test_outer_physical_end_corrections_and_no_inner_boundary_correction():
    stacks, corrections = STAGE08.load_zuber_parameters(RESOURCE)
    assert STAGE08.terminal_dg("AAAA", 4, stacks, corrections) == pytest.approx(3 * -0.94 + 0.22)
    assert STAGE08.terminal_dg("AGAA", 4, stacks, corrections) == pytest.approx(-2.01 - 2.42 - 0.94 + 0.44)
    # Outer G/C terminal pair gets no correction even though the inner sequence contains A/U.
    assert STAGE08.terminal_dg("GAAA", 4, stacks, corrections) == pytest.approx(-2.42 - 0.94 - 0.94)


def test_guide_and_passenger_physical_5p_ends_and_positive_asymmetry_sign():
    stacks, corrections = STAGE08.load_zuber_parameters(RESOURCE)
    passenger = "GCGCAAAA"
    guide = STAGE08.reverse_complement_rna(passenger)
    guide_dg, passenger_dg, ddg = STAGE08.calculate_asymmetry(
        guide, passenger, 4, stacks, corrections
    )
    assert guide.startswith("UUUU") and passenger.startswith("GCGC")
    assert guide_dg == pytest.approx(STAGE08.terminal_dg("UUUU", 4, stacks, corrections))
    assert passenger_dg == pytest.approx(STAGE08.terminal_dg("GCGC", 4, stacks, corrections))
    assert ddg == pytest.approx(guide_dg - passenger_dg)
    assert ddg > 0


def test_no_dicer_overhang_manipulation_is_permitted():
    stacks, corrections = STAGE08.load_zuber_parameters(RESOURCE)
    passenger = "GCGCAAAA"
    guide = STAGE08.reverse_complement_rna(passenger)
    expected = STAGE08.calculate_asymmetry(guide, passenger, 5, stacks, corrections)
    assert expected[0] == STAGE08.terminal_dg(guide, 5, stacks, corrections)
    with pytest.raises(STAGE08.Stage08Error, match="perfect full-length"):
        STAGE08.calculate_asymmetry(guide[2:], passenger, 5, stacks, corrections)


def test_rnafold_parser_preserves_exact_guide_and_finite_structure():
    guide = "AUGCA"
    parsed = STAGE08.parse_rnafold_output(
        ">candidate_a\nAUGCA\n..... ( 0.00)\n", {"candidate_a": guide}
    )
    assert parsed == {"candidate_a": (0.0, ".....")}
    assert len(parsed["candidate_a"][1]) == len(guide)


def test_rnafold_receives_exact_guide_input(monkeypatch):
    captured = {}

    def fake_run(command, input, text, capture_output, check):
        captured["command"] = command
        captured["input"] = input
        return SimpleNamespace(returncode=0, stdout=">c1\nAUGCA\n..... ( 0.00)\n", stderr="")

    monkeypatch.setattr(STAGE08.subprocess, "run", fake_run)
    assert STAGE08.run_rnafold({"c1": "AUGCA"}, 37.0) == {"c1": (0.0, ".....")}
    assert captured["input"] == ">c1\nAUGCA\n"


def test_stage08_has_no_stage07_dependency_or_scoring_columns():
    rule_text = RULE.read_text(encoding="utf-8")
    assert "results/07_" not in rule_text
    forbidden = {"score", "rank", "weight", "threshold", "pass", "fail", "gate"}
    assert not any(
        token in column.lower().split("_")
        for column in STAGE08.BIOPHYSICS_COLUMNS for token in forbidden
    )
    assert len(STAGE08.BIOPHYSICS_COLUMNS) == len(set(STAGE08.BIOPHYSICS_COLUMNS))
