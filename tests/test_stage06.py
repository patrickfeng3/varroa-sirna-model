from __future__ import annotations

import csv
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/stage06.py"
SPEC = importlib.util.spec_from_file_location("stage06", SCRIPT)
assert SPEC and SPEC.loader
stage06 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage06)


def dna_hash(sequence: str) -> str:
    normalized = sequence.upper().replace("U", "T")
    return hashlib.sha256(normalized.encode("ascii")).hexdigest()


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f">{record_id}\n{sequence}\n" for record_id, sequence in records),
        encoding="utf-8",
    )


def manifest_row(
    target_id: str,
    transcript_id: str,
    fasta_name: str,
    record_id: str,
    sequence: str,
    lengths: str,
    annotation: str = "NA",
) -> dict[str, str]:
    normalized = sequence.upper().replace("U", "T")
    return {
        "target_id": target_id,
        "transcript_id": transcript_id,
        "display_name": f"Display {target_id}",
        "organism": "Synthetic organism",
        "molecule_type": "mRNA",
        "fasta_path": fasta_name,
        "fasta_record_id": record_id,
        "annotation_path": annotation,
        "expected_length_nt": str(len(normalized)),
        "sequence_sha256_uppercase_dna": dna_hash(sequence),
        "candidate_lengths_nt": lengths,
        "source_database": "SyntheticDB",
        "source_accession_version": transcript_id,
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=stage06.TARGET_MANIFEST_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize(
    ("input_sequence", "expected"),
    [("acgtt", "ACGTT"), ("ACGUU", "ACGTT"), ("acguu", "ACGTT")],
)
def test_dna_rna_and_lowercase_normalization(input_sequence: str, expected: str) -> None:
    assert stage06.normalize_sequence(input_sequence) == expected


def test_multirecord_fasta_selects_exact_record(tmp_path: Path) -> None:
    fasta = tmp_path / "multi.fa"
    write_fasta(fasta, [("other", "AAAA"), ("wanted", "ccugu"), ("last", "TTTT")])
    assert stage06.load_transcript_sequence(fasta, "wanted") == "CCTGT"
    with pytest.raises(stage06.Stage06Error, match="not found"):
        stage06.load_transcript_sequence(fasta, "missing")


def test_generic_two_target_multiple_length_enumeration_and_global_ids(tmp_path: Path) -> None:
    first_sequence = "acgtuacg"
    second_sequence = "TTGCAAA"
    write_fasta(tmp_path / "first.fa", [("tx_arbitrary", first_sequence)])
    write_fasta(tmp_path / "second.fa", [("tx_second", second_sequence)])
    manifest = tmp_path / "targets.tsv"
    write_manifest(
        manifest,
        [
            manifest_row(
                "toy_alpha", "tx_arbitrary", "first.fa", "tx_arbitrary", first_sequence, "5,6"
            ),
            manifest_row(
                "toy_beta", "tx_second", "second.fa", "tx_second", second_sequence, "4"
            ),
        ],
    )
    output = tmp_path / "results"
    summary = stage06.run_stage06(manifest, output)
    candidates = read_tsv(output / "target_candidates.tsv")

    assert summary["targets"] == 2
    assert summary["candidates"] == (8 - 5 + 1) + (8 - 6 + 1) + (7 - 4 + 1)
    assert len({row["candidate_id"] for row in candidates}) == len(candidates)
    first = next(row for row in candidates if row["candidate_id"] == "toy_alpha__5nt__0001_0005")
    assert first["target_sequence_dna"] == "ACGTT"
    assert first["target_sequence_rna"] == "ACGUU"
    assert first["antisense_guide_sequence_rna"] == "AACGU"
    assert first["annotation_status"] == "unavailable"
    assert first["start_region"] == "NA"


def test_exact_count_first_final_intervals_and_slicing() -> None:
    sequence = "AACCGGTTA"
    target = {"target_id": "generic", "transcript_id": "tx", "display_name": "g", "organism": "o"}
    rows = stage06.enumerate_candidates(target, sequence, [5], None, "unavailable")
    assert len(rows) == len(sequence) - 5 + 1
    assert (rows[0]["start_1based"], rows[0]["end_1based"]) == (1, 5)
    assert (rows[-1]["start_1based"], rows[-1]["end_1based"]) == (5, 9)
    assert rows[0]["target_sequence_dna"] == sequence[:5]
    assert rows[-1]["target_sequence_dna"] == sequence[-5:]


def test_exact_antisense_reverse_complement_and_terminal_orientation() -> None:
    target_rna = "ACGUU"
    guide = stage06.reverse_complement_rna(target_rna)
    assert guide == "AACGU"
    complements = str.maketrans("ACGU", "UGCA")
    assert guide[0] == target_rna[-1].translate(complements)
    assert guide[-1] == target_rna[0].translate(complements)


def test_length_validation_rejects_mismatch_and_oversized_candidate() -> None:
    sequence = "ACGTAC"
    target = manifest_row("x", "tx", "x.fa", "tx", sequence, "5")
    target["expected_length_nt"] = "7"
    with pytest.raises(stage06.Stage06Error, match="length mismatch"):
        stage06.validate_target(target, sequence)
    target["expected_length_nt"] = "6"
    target["candidate_lengths_nt"] = "7"
    with pytest.raises(stage06.Stage06Error, match="exceeds transcript length"):
        stage06.validate_target(target, sequence)


def test_sha256_validation_rejects_mismatch() -> None:
    sequence = "ACGTAC"
    target = manifest_row("x", "tx", "x.fa", "tx", sequence, "5")
    target["sequence_sha256_uppercase_dna"] = "0" * 64
    with pytest.raises(stage06.Stage06Error, match="SHA-256 mismatch"):
        stage06.validate_target(target, sequence)


def test_ambiguous_base_rejected() -> None:
    with pytest.raises(stage06.Stage06Error, match="ambiguous/unexpected"):
        stage06.normalize_sequence("ACNT")


def write_annotation(path: Path, rows: list[tuple[str, str, int, int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(stage06.ANNOTATION_COLUMNS)
        for transcript_id, label, start, end in rows:
            writer.writerow([transcript_id, label, start, end, stage06.COORDINATE_SYSTEM])


def test_complete_annotation_with_arbitrary_labels() -> None:
    regions = [
        {"region_label": "leader", "start_1based": 1, "end_1based": 3},
        {"region_label": "custom_body", "start_1based": 4, "end_1based": 8},
    ]
    annotated = stage06.annotate_candidate(2, 5, regions, "complete")
    assert annotated == {
        "annotation_status": "complete",
        "start_region": "leader",
        "end_region": "custom_body",
        "overlap_regions": "leader;custom_body",
        "crosses_annotation_boundary": "TRUE",
    }


def test_partial_annotation_gap_becomes_unannotated(tmp_path: Path) -> None:
    annotation = tmp_path / "regions.tsv"
    write_annotation(annotation, [("tx", "alpha", 1, 2), ("tx", "omega", 5, 8)])
    regions, status = stage06.load_transcript_regions(annotation, "tx", 8)
    assert status == "partial"
    annotated = stage06.annotate_candidate(2, 6, regions, status)
    assert annotated["overlap_regions"] == "alpha;unannotated;omega"
    assert annotated["start_region"] == "alpha"
    assert annotated["end_region"] == "omega"
    assert annotated["crosses_annotation_boundary"] == "TRUE"


def test_annotation_unavailable_fields_are_na() -> None:
    regions, status = stage06.load_transcript_regions(None, "anything", 12)
    assert regions is None
    assert status == "unavailable"
    annotated = stage06.annotate_candidate(1, 5, regions, status)
    assert annotated == {
        "annotation_status": "unavailable",
        "start_region": "NA",
        "end_region": "NA",
        "overlap_regions": "NA",
        "crosses_annotation_boundary": "NA",
    }


def test_overlapping_annotation_is_rejected(tmp_path: Path) -> None:
    annotation = tmp_path / "regions.tsv"
    write_annotation(annotation, [("tx", "left", 1, 5), ("tx", "right", 5, 8)])
    with pytest.raises(stage06.Stage06Error, match="overlap"):
        stage06.load_transcript_regions(annotation, "tx", 8)


def test_vd_chibin_reference_identity_counts_boundaries_and_orientation(tmp_path: Path) -> None:
    manifest = ROOT / "resources/targets/target_manifest.tsv"
    output = tmp_path / "stage06"
    stage06.run_stage06(manifest, output)
    references = read_tsv(output / "target_reference_summary.tsv")
    candidates = read_tsv(output / "target_candidates.tsv")

    assert len(references) == 1
    assert references[0]["normalized_length_nt"] == "710"
    assert references[0]["sequence_sha256_uppercase_dna"] == (
        "4a0d25aa05b269a118ed1b952dca63ccd1c0a7978fc42295faf3bf650e43ea42"
    )
    by_length = {
        length: [row for row in candidates if row["candidate_length_nt"] == str(length)]
        for length in (23, 24)
    }
    assert len(by_length[23]) == 688
    assert len(by_length[24]) == 687
    assert len(candidates) == 1375

    for length, expected_boundary_count in ((23, 22), (24, 23)):
        rows = by_length[length]
        assert sum(
            row["start_region"] == "5_prime_UTR" and row["end_region"] == "CDS"
            for row in rows
        ) == expected_boundary_count
        assert sum(
            row["start_region"] == "CDS" and row["end_region"] == "3_prime_UTR"
            for row in rows
        ) == expected_boundary_count

    complement = str.maketrans("ACGU", "UGCA")
    for row in (candidates[0], candidates[len(candidates) // 2], candidates[-1]):
        target_rna = row["target_sequence_rna"]
        guide = row["antisense_guide_sequence_rna"]
        assert guide == stage06.reverse_complement_rna(target_rna)
        assert guide[0] == target_rna[-1].translate(complement)
        assert guide[-1] == target_rna[0].translate(complement)


def test_cli_writes_exact_four_outputs(tmp_path: Path) -> None:
    sequence = "ACGTTACG"
    write_fasta(tmp_path / "toy.fa", [("toy_tx", sequence)])
    manifest = tmp_path / "targets.tsv"
    write_manifest(
        manifest,
        [manifest_row("toy", "toy_tx", "toy.fa", "toy_tx", sequence, "5")],
    )
    output = tmp_path / "output"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--target-manifest",
            str(manifest),
            "--output-root",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout
    assert sorted(path.relative_to(output).as_posix() for path in output.rglob("*.tsv")) == [
        "provenance/stage06_manifest.tsv",
        "qc/stage06_accounting.tsv",
        "target_candidates.tsv",
        "target_reference_summary.tsv",
    ]
