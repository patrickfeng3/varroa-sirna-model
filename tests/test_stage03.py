import importlib.util
import shutil
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "workflow/scripts/stage03.py"
SPEC = importlib.util.spec_from_file_location("stage03", SCRIPT)
stage03 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = stage03
SPEC.loader.exec_module(stage03)


CONFIG = stage03.Stage03Config((23, 24), 15, 30, 2, -2, "1.0.6")


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


def find_run(result, length=23, strand="sense"):
    return next(
        row for row in result["runs"]
        if row["sample"] == "S1" and row["analysis_unit"] == "V1"
        and row["focal_length"] == length and row["focal_strand"] == strand
    )


def test_file_a_collapses_duplicate_sequence_and_sums_focal_abundance():
    rows = [feature(count=7), feature(count=13)]
    result = stage03.collapse_inputs([eligibility()], rows, CONFIG)
    run = find_run(result)
    focals = result["focal_by_run"][run["run_id"]]
    assert len(focals) == 1
    assert focals[0]["focal_abundance"] == 20


def test_file_b_collapses_and_uses_15_30_inclusive_range():
    rows = [
        feature(),
        feature(length=15, strand="antisense", sequence="C" * 15, count=2),
        feature(length=15, strand="antisense", sequence="C" * 15, count=3),
        feature(length=30, strand="antisense", sequence="G" * 30),
        feature(length=14, strand="antisense", sequence="T" * 14),
        feature(length=31, strand="antisense", sequence="A" * 31),
    ]
    result = stage03.collapse_inputs([eligibility()], rows, CONFIG)
    run = find_run(result)
    passengers = result["passenger_by_run"][run["run_id"]]
    assert sorted(row["passenger_length"] for row in passengers) == [15, 30]
    assert next(row for row in passengers if row["passenger_length"] == 15)["passenger_abundance"] == 5


def test_file_b_is_same_sample_virus_and_opposite_mapped_strand():
    rows = [
        feature(),
        feature(strand="antisense", sequence="C" * 23),
        feature(strand="sense", sequence="G" * 23),
        feature(unit="V2", strand="antisense", sequence="T" * 23),
        feature(sample="S2", strand="antisense", sequence="AC" * 11 + "A"),
    ]
    result = stage03.collapse_inputs([eligibility(), eligibility("S1", "V2"), eligibility("S2", "V1")], rows, CONFIG)
    run = find_run(result)
    assert [row["sequence"] for row in result["passenger_by_run"][run["run_id"]]] == ["C" * 23]


def test_file_a_and_b_preserve_observed_physical_orientation():
    focal = "ACGTTGCACTGATCGTACGATGC"
    passenger = "ATCGTACGATCAGTGCAACGTAA"
    result = stage03.collapse_inputs(
        [eligibility()],
        [feature(sequence=focal), feature(strand="antisense", sequence=passenger)],
        CONFIG,
    )
    run = find_run(result)
    assert result["focal_by_run"][run["run_id"]][0]["sequence"] == focal
    assert result["passenger_by_run"][run["run_id"]][0]["sequence"] == passenger


def test_stable_fasta_identifiers_are_order_independent_and_unique():
    rows = [feature(sequence="A" * 23), feature(sequence="C" * 23)]
    first = stage03.collapse_inputs([eligibility()], rows, CONFIG)
    second = stage03.collapse_inputs([eligibility()], reversed(rows), CONFIG)
    ids1 = sorted(row["focal_id"] for row in first["focal_manifest"])
    ids2 = sorted(row["focal_id"] for row in second["focal_manifest"])
    assert ids1 == ids2
    assert len(ids1) == len(set(ids1))


@pytest.fixture(scope="session")
def official_preflight(tmp_path_factory):
    root = tmp_path_factory.mktemp("official_preflight")
    result = stage03.run_preflight(root, "1.0.6")
    return root, result


def synthetic_run_and_manifests():
    run = {
        "run_id": "synthetic", "sample": "S1", "analysis_unit": "V1",
        "biological_virus": "V1", "focal_length": 23, "focal_strand": "sense",
    }
    refs, reads = stage03.synthetic_records()
    abundances = {"FJ": 10, "FB": 20, "FO": 30}
    focals = [
        {"focal_id": identifier, "sequence": sequence, "focal_abundance": abundances[identifier]}
        for identifier, sequence in refs
    ]
    passengers = [
        {"passenger_id": identifier, "sequence": sequence}
        for identifier, sequence in reads
    ]
    return run, focals, passengers


def test_official_synthetic_preflight_produces_required_outputs(official_preflight):
    root, result = official_preflight
    assert result["versions"]["stepRNA"] == "1.0.6"
    assert all(row["status"] == "PASS" for row in result["checks"])
    assert (root / "provenance/preflight/official/synthetic_overhang.csv").exists()
    assert (root / "provenance/preflight/official/synthetic.sorted.bam").exists()


def test_official_sign_convention_and_known_joint_geometry(official_preflight):
    _, result = official_preflight
    assert result["geometries"][("FO", "QO")][0] == -2
    assert result["geometries"][("FB", "QB")] == (0, 0)
    assert result["geometries"][("FJ", "QJ")] == (2, -2)


@pytest.fixture
def synthetic_parsed(official_preflight):
    root, _ = official_preflight
    run, focals, passengers = synthetic_run_and_manifests()
    raw = root / "provenance/preflight/official"
    return stage03.parse_official_run(run, raw, focals, passengers, CONFIG)


def test_parser_keeps_official_counts_separate_from_focal_abundance(synthetic_parsed):
    five_prime = [row for row in synthetic_parsed["spectrum"] if row["end"] == "5p"]
    assert sum(row["official_duplex_count"] for row in five_prime) == 3
    assert synthetic_parsed["recovery"]["total_focal_abundance"] == 60


def test_passenger_recovery_fraction_unique(synthetic_parsed):
    assert synthetic_parsed["recovery"]["passenger_recovery_fraction_unique"] == 1


def test_passenger_recovery_fraction_abundance(synthetic_parsed):
    assert synthetic_parsed["recovery"]["passenger_recovery_fraction_abundance"] == 1


def test_zero_passenger_recovery_behavior():
    run, focals, _ = synthetic_run_and_manifests()
    recovery, geometry, _ = stage03.summarise_duplexes(run, focals, [], CONFIG)
    assert recovery["passenger_recovery_fraction_unique"] == 0
    assert recovery["passenger_recovery_fraction_abundance"] == 0
    assert geometry["varroa_2nt_joint_duplex_fraction"] is None
    assert geometry["varroa_2nt_reference_fraction_recovered"] is None


def geometry_fixture():
    run, focals, _ = synthetic_run_and_manifests()
    duplexes = [
        {"focal_id": "FJ", "steprna_5p_distance": 2, "steprna_3p_distance": 0},
        {"focal_id": "FJ", "steprna_5p_distance": 0, "steprna_3p_distance": -2},
        {"focal_id": "FB", "steprna_5p_distance": 2, "steprna_3p_distance": -2},
        {"focal_id": "FB", "steprna_5p_distance": -1, "steprna_3p_distance": 1},
    ]
    return stage03.summarise_duplexes(run, focals, duplexes, CONFIG)


def test_joint_geometry_requires_both_distances_from_same_duplex():
    _, geometry, _ = geometry_fixture()
    assert geometry["n_joint_geometry_duplexes"] == 1


def test_joint_duplex_fraction():
    _, geometry, _ = geometry_fixture()
    assert geometry["varroa_2nt_joint_duplex_fraction"] == 0.25


def test_joint_reference_fraction_all_focals():
    _, geometry, _ = geometry_fixture()
    assert geometry["varroa_2nt_reference_fraction_all"] == 1 / 3


def test_joint_reference_fraction_among_recovered_focals():
    _, geometry, _ = geometry_fixture()
    assert geometry["varroa_2nt_reference_fraction_recovered"] == 0.5


def test_abundance_weighted_joint_reference_support_fractions():
    _, geometry, _ = geometry_fixture()
    assert geometry["varroa_2nt_reference_fraction_abundance_all"] == 20 / 60
    assert geometry["varroa_2nt_reference_fraction_abundance_recovered"] == 20 / 30


def test_one_focal_reference_may_support_multiple_geometries():
    _, geometry, refs = geometry_fixture()
    assert geometry["n_focal_references_supporting_joint_geometry"] == 1
    assert refs[0]["n_recovered_duplexes_for_reference"] == 2
    assert refs[0]["n_joint_geometry_duplexes_for_reference"] == 1


def test_passenger_length_parsing(synthetic_parsed):
    counts = {
        row["passenger_length"]: row["official_duplex_count"]
        for row in synthetic_parsed["passenger_lengths"]
    }
    assert counts == {21: 1, 23: 2}


def test_missing_official_output_causes_structured_failure(tmp_path):
    run, focals, passengers = synthetic_run_and_manifests()
    with pytest.raises(stage03.Stage03Error, match="missing official files"):
        stage03.parse_official_run(run, tmp_path, focals, passengers, CONFIG)


def test_malformed_official_output_causes_structured_failure(official_preflight, tmp_path):
    root, _ = official_preflight
    raw = tmp_path / "official"
    shutil.copytree(root / "provenance/preflight/official", raw)
    (raw / "synthetic_overhang.csv").write_text("bad,header\n1,2\n")
    run, focals, passengers = synthetic_run_and_manifests()
    with pytest.raises(stage03.Stage03Error, match="malformed official overhang"):
        stage03.parse_official_run(run, raw, focals, passengers, CONFIG)


def test_parser_ids_must_exist_in_manifests(official_preflight):
    root, _ = official_preflight
    run, focals, passengers = synthetic_run_and_manifests()
    with pytest.raises(stage03.Stage03Error, match="identifier mismatch"):
        stage03.parse_official_run(
            run, root / "provenance/preflight/official", focals[1:], passengers, CONFIG
        )
