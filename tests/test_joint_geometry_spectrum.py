import importlib.util
import sys
from pathlib import Path


SCRIPTS = Path(__file__).parents[1] / "workflow/scripts"


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


stage03 = load("stage03", "stage03.py")
stage03_joint = load("stage03_joint_spectrum", "stage03_joint_spectrum.py")
stage04 = load("stage04", "stage04.py")
stage04_joint = load("stage04_joint_spectrum", "stage04_joint_spectrum.py")

CONFIG = stage04.Stage04Config(100, 20260814, "percentile", 0.95, 1e-12, 2, -2)


def duplex(d5, d3):
    return {"steprna_5p_distance": d5, "steprna_3p_distance": d3}


def pair_row(d5, d3, count, total=4, run_id="R1"):
    return {
        "sample": "S1", "analysis_unit": "V1", "biological_virus": "V1",
        "focal_length": 23, "focal_strand": "antisense", "run_id": run_id,
        "steprna_5p_distance": d5, "steprna_3p_distance": d3,
        "official_duplex_count": count, "total_recovered_duplexes": total,
        "joint_duplex_fraction": count / total,
    }


def test_same_duplex_joint_counting_does_not_combine_marginals():
    rows = stage03_joint.joint_geometry_spectrum([
        duplex(0, -2), duplex(2, 0), duplex(0, 0), duplex(2, -2),
    ])
    counts = {
        (row["steprna_5p_distance"], row["steprna_3p_distance"]): row["official_duplex_count"]
        for row in rows
    }
    assert counts == {(0, -2): 1, (0, 0): 1, (2, -2): 1, (2, 0): 1}


def test_nonempty_joint_fractions_sum_to_one():
    rows = stage03_joint.joint_geometry_spectrum([
        duplex(0, 0), duplex(0, 0), duplex(2, -2), duplex(-1, 3),
    ])
    assert sum(row["official_duplex_count"] for row in rows) == 4
    assert abs(sum(row["joint_duplex_fraction"] for row in rows) - 1.0) < 1e-12


def test_known_zero_zero_and_plus2_minus2_geometries_are_retained_and_aggregated():
    sparse = [pair_row(0, 0, 2), pair_row(2, -2, 1), pair_row(-1, 3, 1)]
    _, _, across, modes, summary, qc = stage04_joint.aggregate_joint_spectrum(sparse, CONFIG)
    lookup = {
        (row["steprna_5p_distance"], row["steprna_3p_distance"]): row
        for row in across
    }
    assert lookup[(0, 0)]["sample_balanced_median_joint_duplex_fraction"] == 0.5
    assert lookup[(2, -2)]["sample_balanced_median_joint_duplex_fraction"] == 0.25
    assert modes[0]["zero_zero_is_most_common"] == 1
    assert summary[0]["n_runs_zero_zero_fraction_gt_0_5"] == 0
    assert not [row for row in qc if row["status"] == "FAIL"]


def test_prespecified_plus2_minus2_output_regression_guard():
    spectrum = stage03_joint.joint_geometry_spectrum([
        duplex(0, 0), duplex(2, -2), duplex(2, -2), duplex(1, -1),
    ])
    assert stage03_joint.agrees_with_prespecified(spectrum, 4, 2)
    assert not stage03_joint.agrees_with_prespecified(spectrum, 4, 1)
    assert not stage03_joint.agrees_with_prespecified(spectrum, 5, 2)
