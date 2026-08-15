#!/usr/bin/env python3
"""Build the Stage 03 full same-duplex geometry spectrum from existing BAMs."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Iterable

import pysam

from stage03 import (
    official_distances,
    write_table,
)


SPECTRUM_FIELDS = [
    "sample", "analysis_unit", "biological_virus", "focal_length",
    "focal_strand", "steprna_5p_distance", "steprna_3p_distance",
    "official_duplex_count", "total_recovered_duplexes",
    "joint_duplex_fraction", "run_id",
]
QC_FIELDS = ["metric", "status", "value", "details"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def joint_geometry_spectrum(
    duplexes: Iterable[dict[str, object]],
) -> list[dict[str, object]]:
    """Count complete same-duplex 5p/3p geometries without combining marginals."""
    counts = Counter(
        (int(row["steprna_5p_distance"]), int(row["steprna_3p_distance"]))
        for row in duplexes
    )
    total = sum(counts.values())
    return [
        {
            "steprna_5p_distance": d5,
            "steprna_3p_distance": d3,
            "official_duplex_count": count,
            "total_recovered_duplexes": total,
            "joint_duplex_fraction": count / total,
        }
        for (d5, d3), count in sorted(counts.items())
    ]


def agrees_with_prespecified(
    spectrum: list[dict[str, object]], expected_total: int, expected_joint: int
) -> bool:
    """Protect the existing pre-specified (+2,-2) result from parser regression."""
    observed_total = sum(int(row["official_duplex_count"]) for row in spectrum)
    observed_joint = sum(
        int(row["official_duplex_count"])
        for row in spectrum
        if int(row["steprna_5p_distance"]) == 2
        and int(row["steprna_3p_distance"]) == -2
    )
    return observed_total == expected_total and observed_joint == expected_joint


def build_joint_spectrum(stage03_root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    runs = read_tsv(stage03_root / "provenance/run_manifest.tsv")
    existing_joint = {
        row["run_id"]: row
        for row in read_tsv(stage03_root / "parsed/joint_geometry_by_pair.tsv")
    }
    output: list[dict[str, object]] = []
    sum_inconsistencies = 0
    prespecified_regressions = 0
    missing_bams = 0
    parsed_runs = 0
    for run in runs:
        if run["status"] != "success":
            continue
        run_id = run["run_id"]
        passed = (
            stage03_root / "raw" / run_id / f"{run_id}_AlignmentFiles"
            / f"{run_id}_passed.bam"
        )
        if not passed.exists():
            missing_bams += 1
            continue
        duplexes = []
        with pysam.AlignmentFile(passed, "rb") as bam:
            for alignment in bam:
                d5, d3 = official_distances(bam, alignment)
                duplexes.append({
                    "steprna_5p_distance": d5,
                    "steprna_3p_distance": d3,
                })
        spectrum = joint_geometry_spectrum(duplexes)
        total = len(duplexes)
        sum_inconsistencies += int(
            sum(int(row["official_duplex_count"]) for row in spectrum) != total
        )
        expected = existing_joint.get(run_id)
        if (
            expected is None
            or not agrees_with_prespecified(
                spectrum,
                int(expected["n_recovered_duplexes"]),
                int(expected["n_joint_geometry_duplexes"]),
            )
        ):
            prespecified_regressions += 1
        common = {
            "sample": run["sample"],
            "analysis_unit": run["analysis_unit"],
            "biological_virus": run["biological_virus"],
            "focal_length": int(run["focal_length"]),
            "focal_strand": run["focal_strand"],
            "run_id": run_id,
        }
        output.extend({**common, **row} for row in spectrum)
        parsed_runs += 1
    qc = [
        {"metric": "successful_runs_expected", "status": "PASS", "value": sum(r["status"] == "success" for r in runs), "details": ""},
        {"metric": "joint_spectrum_runs_parsed", "status": "PASS", "value": parsed_runs, "details": ""},
        {"metric": "missing_passed_bams", "status": "FAIL" if missing_bams else "PASS", "value": missing_bams, "details": ""},
        {"metric": "joint_geometry_count_sum_inconsistencies", "status": "FAIL" if sum_inconsistencies else "PASS", "value": sum_inconsistencies, "details": ""},
        {"metric": "prespecified_plus2_minus2_regressions", "status": "FAIL" if prespecified_regressions else "PASS", "value": prespecified_regressions, "details": ""},
    ]
    return output, qc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage03-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    args = parser.parse_args(argv)
    started = time.monotonic()
    rows, qc = build_joint_spectrum(args.stage03_root.resolve())
    write_table(args.output.resolve(), rows, SPECTRUM_FIELDS)
    write_table(args.qc_output.resolve(), qc, QC_FIELDS)
    print(f"Stage 03 joint-spectrum post-processing completed in {time.monotonic()-started:.3f} seconds", file=sys.stderr)
    return 1 if any(row["status"] == "FAIL" for row in qc) else 0


if __name__ == "__main__":
    sys.exit(main())
