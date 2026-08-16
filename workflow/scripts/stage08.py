#!/usr/bin/env python3
"""Canonical Stage 08: generic raw candidate-biophysics descriptors."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Iterable, Sequence


STAGE06_COLUMNS = [
    "target_id", "transcript_id", "display_name", "organism", "candidate_id",
    "candidate_length_nt", "start_1based", "end_1based", "target_sequence_dna",
    "target_sequence_rna", "antisense_guide_sequence_rna", "annotation_status",
    "start_region", "end_region", "overlap_regions", "crosses_annotation_boundary",
]

BIOPHYSICS_COLUMNS = STAGE06_COLUMNS + [
    "target_whole_p_unpaired", "target_whole_p_unpaired_w100_l80",
    "target_whole_p_unpaired_w200_l150", "target_seed_g2_8_p_unpaired",
    "target_seed_g2_8_p_unpaired_w100_l80",
    "target_seed_g2_8_p_unpaired_w200_l150", "guide_5p_terminal_dg_4bp",
    "passenger_5p_terminal_dg_4bp", "asymmetry_ddg_4bp",
    "guide_5p_terminal_dg_5bp", "passenger_5p_terminal_dg_5bp",
    "asymmetry_ddg_5bp", "guide_self_fold_mfe_kcal_mol",
    "guide_self_fold_structure",
]

QC_COLUMNS = ["status", "check", "target_id", "candidate_length_nt", "value", "details"]
PARAMETER_COLUMNS = [
    "target_id", "transcript_id", "parameter_set", "viennarna_version",
    "temperature_c", "rnaplfold_window_nt", "rnaplfold_max_bp_span_nt",
    "rnaplfold_ulength_nt", "rnaplfold_effective_window_nt",
    "rnaplfold_effective_max_bp_span_nt", "terminal_thermo_parameter_source",
    "terminal_thermo_parameter_resource", "terminal_thermo_parameter_resource_sha256",
    "candidate_duplex_overhang_assumption", "self_fold_method",
    "stage06_candidate_table_sha256", "target_manifest_sha256",
    "transcript_sequence_sha256_uppercase_dna",
]

MANIFEST_REQUIRED = {
    "target_id", "transcript_id", "fasta_path", "fasta_record_id",
    "expected_length_nt", "sequence_sha256_uppercase_dna",
}
ZUBER_COLUMNS = [
    "parameter_type", "parameter_id", "source_stack", "delta_g37_kcal_mol",
    "units", "source_doi", "description",
]
RNA_ALPHABET = set("ACGU")
DNA_ALPHABET = set("ACGT")
ZUBER_DOI = "10.1093/nar/gkac261"
ACCESSIBILITY_OUTPUTS = {
    "W150_L100_main": ("target_whole_p_unpaired", "target_seed_g2_8_p_unpaired"),
    "W100_L80_sensitivity": (
        "target_whole_p_unpaired_w100_l80", "target_seed_g2_8_p_unpaired_w100_l80"
    ),
    "W200_L150_sensitivity": (
        "target_whole_p_unpaired_w200_l150", "target_seed_g2_8_p_unpaired_w200_l150"
    ),
}


class Stage08Error(ValueError):
    """Deterministic Stage 08 validation or execution error."""


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return reader.fieldnames or [], list(reader)
    except OSError as exc:
        raise Stage08Error(f"cannot read TSV {path}: {exc}") from exc


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _format_value(row.get(field)) for field in fields})
    os.replace(temporary, path)


def _format_value(value: object) -> object:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return value


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_sha256_dna(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def normalize_dna(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    if not normalized or set(normalized) - DNA_ALPHABET:
        raise Stage08Error("transcript sequence contains empty or ambiguous/non-ACGT bases")
    return normalized


def normalize_rna(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper().replace("T", "U")
    if not normalized or set(normalized) - RNA_ALPHABET:
        raise Stage08Error("RNA sequence contains empty or ambiguous/non-ACGU bases")
    return normalized


def reverse_complement_rna(sequence: str) -> str:
    return normalize_rna(sequence).translate(str.maketrans("ACGU", "UGCA"))[::-1]


def read_fasta_record(path: Path, record_id: str) -> str:
    records: dict[str, list[str]] = {}
    current: str | None = None
    try:
        with path.open(encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    if current in records:
                        raise Stage08Error(f"duplicate FASTA record ID {current} in {path}")
                    records[current] = []
                elif current is None:
                    raise Stage08Error(f"sequence before FASTA header in {path}")
                else:
                    records[current].append(line)
    except OSError as exc:
        raise Stage08Error(f"cannot read FASTA {path}: {exc}") from exc
    if record_id not in records:
        raise Stage08Error(f"FASTA record {record_id} not found in {path}")
    return normalize_dna("".join(records[record_id]))


def _resolve_manifest_path(manifest_path: Path, value: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    repository_relative = Path.cwd() / candidate
    return repository_relative.resolve() if repository_relative.exists() else (manifest_path.parent / candidate).resolve()


def load_registered_transcripts(manifest_path: Path) -> dict[tuple[str, str], dict[str, object]]:
    header, rows = _read_tsv(manifest_path)
    missing = MANIFEST_REQUIRED - set(header)
    if missing:
        raise Stage08Error(f"target manifest missing columns: {sorted(missing)}")
    transcripts: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        target_id = row["target_id"].strip()
        transcript_id = row["transcript_id"].strip()
        key = (target_id, transcript_id)
        if key in transcripts:
            raise Stage08Error(f"duplicate target/transcript manifest row: {key}")
        fasta_path = _resolve_manifest_path(manifest_path, row["fasta_path"].strip())
        sequence_dna = read_fasta_record(fasta_path, row["fasta_record_id"].strip())
        expected_length = int(row["expected_length_nt"])
        expected_sha = row["sequence_sha256_uppercase_dna"].strip().lower()
        observed_sha = sequence_sha256_dna(sequence_dna)
        if len(sequence_dna) != expected_length:
            raise Stage08Error(f"registered length mismatch for {key}: {len(sequence_dna)} != {expected_length}")
        if observed_sha != expected_sha:
            raise Stage08Error(f"registered SHA-256 mismatch for {key}: {observed_sha} != {expected_sha}")
        transcripts[key] = {
            "target_id": target_id, "transcript_id": transcript_id,
            "sequence_dna": sequence_dna, "sequence_rna": sequence_dna.replace("T", "U"),
            "sequence_sha256": observed_sha, "fasta_path": str(fasta_path),
        }
    if not transcripts:
        raise Stage08Error("target manifest contains no registered transcripts")
    return transcripts


def load_and_validate_candidates(
    candidate_path: Path, transcripts: dict[tuple[str, str], dict[str, object]]
) -> list[dict[str, object]]:
    header, source_rows = _read_tsv(candidate_path)
    if header != STAGE06_COLUMNS:
        raise Stage08Error(f"Stage 06 candidate schema mismatch: observed {header}")
    if not source_rows:
        raise Stage08Error("Stage 06 candidate table contains no candidates")
    seen_ids: set[str] = set()
    candidates: list[dict[str, object]] = []
    for line_number, source in enumerate(source_rows, start=2):
        row: dict[str, object] = dict(source)
        key = (source["target_id"], source["transcript_id"])
        if key not in transcripts:
            raise Stage08Error(f"candidate line {line_number} has unregistered target/transcript {key}")
        candidate_id = source["candidate_id"]
        if candidate_id in seen_ids:
            raise Stage08Error(f"duplicate candidate_id {candidate_id}")
        seen_ids.add(candidate_id)
        length = int(source["candidate_length_nt"])
        start = int(source["start_1based"])
        end = int(source["end_1based"])
        transcript_dna = str(transcripts[key]["sequence_dna"])
        if length <= 0 or end - start + 1 != length or not (1 <= start <= end <= len(transcript_dna)):
            raise Stage08Error(f"invalid candidate coordinates/length at line {line_number}")
        target_dna = normalize_dna(source["target_sequence_dna"])
        target_rna = normalize_rna(source["target_sequence_rna"])
        guide_rna = normalize_rna(source["antisense_guide_sequence_rna"])
        expected_dna = transcript_dna[start - 1:end]
        if target_dna != expected_dna or target_rna != expected_dna.replace("T", "U"):
            raise Stage08Error(f"candidate transcript slice mismatch for {candidate_id}")
        if len(target_rna) != length or len(guide_rna) != length:
            raise Stage08Error(f"candidate sequence length mismatch for {candidate_id}")
        if guide_rna != reverse_complement_rna(target_rna):
            raise Stage08Error(f"candidate guide reverse-complement mismatch for {candidate_id}")
        row.update({
            "candidate_length_nt": length, "start_1based": start, "end_1based": end,
            "target_sequence_dna": target_dna, "target_sequence_rna": target_rna,
            "antisense_guide_sequence_rna": guide_rna,
        })
        candidates.append(row)
    return candidates


class LunpTable:
    """RNAplfold unpaired probabilities indexed as (interval end, interval length)."""

    def __init__(self, probabilities: dict[tuple[int, int], float]):
        self.probabilities = probabilities

    def interval_probability(self, start_1based: int, end_1based: int) -> float:
        if start_1based < 1 or end_1based < start_1based:
            raise Stage08Error(f"invalid unpaired interval [{start_1based},{end_1based}]")
        length = end_1based - start_1based + 1
        key = (end_1based, length)
        if key not in self.probabilities:
            raise Stage08Error(
                f"RNAplfold _lunp value absent for interval [{start_1based},{end_1based}] "
                f"(row/end={end_1based}, length-column={length})"
            )
        value = self.probabilities[key]
        validate_probability(value)
        return value


def parse_lunp(path: Path) -> LunpTable:
    probabilities: dict[tuple[int, int], float] = {}
    try:
        with path.open(encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                try:
                    end = int(fields[0])
                except (ValueError, IndexError) as exc:
                    raise Stage08Error(f"malformed RNAplfold _lunp row: {line}") from exc
                for length, text in enumerate(fields[1:], start=1):
                    if text.upper() in {"NA", "NAN"}:
                        continue
                    try:
                        value = float(text)
                    except ValueError as exc:
                        raise Stage08Error(f"malformed RNAplfold probability: {text}") from exc
                    validate_probability(value)
                    probabilities[(end, length)] = value
    except OSError as exc:
        raise Stage08Error(f"cannot read RNAplfold output {path}: {exc}") from exc
    if not probabilities:
        raise Stage08Error(f"RNAplfold output contains no probabilities: {path}")
    return LunpTable(probabilities)


def validate_probability(value: float) -> float:
    if not math.isfinite(value) or value < 0 or value > 1:
        raise Stage08Error(f"invalid accessibility probability: {value}")
    return value


def effective_window_parameters(requested_w: int, requested_l: int, transcript_length: int) -> tuple[int, int]:
    if requested_w <= 0 or requested_l <= 0 or transcript_length <= 1:
        raise Stage08Error("RNAplfold W/L and transcript length must be positive; transcript length must exceed 1")
    effective_w = min(requested_w, transcript_length)
    effective_l = min(requested_l, effective_w - 1)
    return effective_w, effective_l


def run_rnaplfold(
    transcript_rna: str, window_nt: int, max_bp_span_nt: int, ulength_nt: int,
    temperature_c: float, executable: str = "RNAplfold",
) -> LunpTable:
    with tempfile.TemporaryDirectory(prefix="stage08_plfold_") as directory:
        workdir = Path(directory)
        command = [
            executable, "-W", str(window_nt), "-L", str(max_bp_span_nt),
            "-u", str(ulength_nt), "-T", str(temperature_c),
        ]
        completed = subprocess.run(
            command, input=normalize_rna(transcript_rna) + "\n", text=True,
            cwd=workdir, capture_output=True, check=False,
        )
        if completed.returncode != 0:
            raise Stage08Error(f"RNAplfold failed ({completed.returncode}): {completed.stderr.strip()}")
        lunp_files = sorted(set(workdir.glob("*_lunp")) | set(workdir.glob("plfold_lunp")))
        if len(lunp_files) != 1:
            raise Stage08Error(f"RNAplfold produced {len(lunp_files)} _lunp files; expected exactly one")
        return parse_lunp(lunp_files[0])


def calculate_accessibilities(
    candidates: Sequence[dict[str, object]],
    transcripts: dict[tuple[str, str], dict[str, object]],
    parameter_sets: Sequence[dict[str, object]],
    temperature_c: float,
    fold_runner: Callable[[str, int, int, int, float], LunpTable] = run_rnaplfold,
) -> tuple[dict[str, dict[str, float | None]], list[dict[str, object]], int]:
    by_transcript: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for candidate in candidates:
        by_transcript[(str(candidate["target_id"]), str(candidate["transcript_id"]))].append(candidate)
    values: dict[str, dict[str, float | None]] = {str(row["candidate_id"]): {} for row in candidates}
    provenance: list[dict[str, object]] = []
    run_count = 0
    for key in sorted(by_transcript):
        transcript_candidates = by_transcript[key]
        transcript_rna = str(transcripts[key]["sequence_rna"])
        ulength = max(max(int(row["candidate_length_nt"]) for row in transcript_candidates), 7)
        for parameter in parameter_sets:
            parameter_id = str(parameter["id"])
            if parameter_id not in ACCESSIBILITY_OUTPUTS:
                raise Stage08Error(f"unsupported Stage 08 accessibility parameter ID: {parameter_id}")
            requested_w = int(parameter["window_nt"])
            requested_l = int(parameter["max_bp_span_nt"])
            effective_w, effective_l = effective_window_parameters(
                requested_w, requested_l, len(transcript_rna)
            )
            table = fold_runner(transcript_rna, effective_w, effective_l, ulength, temperature_c)
            run_count += 1
            whole_field, seed_field = ACCESSIBILITY_OUTPUTS[parameter_id]
            for candidate in transcript_candidates:
                start = int(candidate["start_1based"])
                end = int(candidate["end_1based"])
                candidate_id = str(candidate["candidate_id"])
                values[candidate_id][whole_field] = table.interval_probability(start, end)
                values[candidate_id][seed_field] = (
                    table.interval_probability(end - 7, end - 1)
                    if int(candidate["candidate_length_nt"]) >= 8 else None
                )
            provenance.append({
                "target_id": key[0], "transcript_id": key[1], "parameter_set": parameter_id,
                "temperature_c": temperature_c, "rnaplfold_window_nt": requested_w,
                "rnaplfold_max_bp_span_nt": requested_l, "rnaplfold_ulength_nt": ulength,
                "rnaplfold_effective_window_nt": effective_w,
                "rnaplfold_effective_max_bp_span_nt": effective_l,
            })
    return values, provenance, run_count


def load_zuber_parameters(path: Path) -> tuple[dict[str, float], dict[str, float]]:
    header, rows = _read_tsv(path)
    if header != ZUBER_COLUMNS:
        raise Stage08Error(f"Zuber resource schema mismatch: observed {header}")
    stacks: dict[str, float] = {}
    corrections: dict[str, float] = {}
    for row in rows:
        if row["units"] != "kcal/mol" or row["source_doi"] != ZUBER_DOI:
            raise Stage08Error("Zuber resource has incorrect units or source DOI")
        value = float(row["delta_g37_kcal_mol"])
        if not math.isfinite(value):
            raise Stage08Error("Zuber resource contains non-finite parameter")
        if row["parameter_type"] == "stack":
            key = row["parameter_id"]
            if len(key) != 2 or set(key) - RNA_ALPHABET or key in stacks:
                raise Stage08Error(f"invalid/duplicate Zuber stack key: {key}")
            stacks[key] = value
        elif row["parameter_type"] == "end_correction":
            key = row["parameter_id"]
            if key in corrections:
                raise Stage08Error(f"duplicate Zuber end-correction key: {key}")
            corrections[key] = value
        else:
            raise Stage08Error(f"unsupported Zuber parameter type: {row['parameter_type']}")
    expected_stacks = {left + right for left in "ACGU" for right in "ACGU"}
    expected_corrections = {
        "AU_terminal_on_AU_penultimate", "AU_terminal_on_GC_penultimate"
    }
    if set(stacks) != expected_stacks or set(corrections) != expected_corrections:
        raise Stage08Error("Zuber resource does not contain exactly 16 WCF stacks and two end corrections")
    for stack, value in stacks.items():
        reverse_complement = reverse_complement_rna(stack)
        if value != stacks[reverse_complement]:
            raise Stage08Error(f"reverse-complement-equivalent Zuber values disagree: {stack}")
    return stacks, corrections


def terminal_dg(
    strand_sequence_rna: str, terminal_length_nt: int,
    stacks: dict[str, float], corrections: dict[str, float],
) -> float:
    sequence = normalize_rna(strand_sequence_rna)
    if terminal_length_nt < 2 or len(sequence) < terminal_length_nt:
        raise Stage08Error("terminal stability window is invalid for strand length")
    segment = sequence[:terminal_length_nt]
    stack_sum = sum(stacks[segment[index:index + 2]] for index in range(terminal_length_nt - 1))
    correction = 0.0
    if segment[0] in "AU":
        correction_key = (
            "AU_terminal_on_AU_penultimate"
            if segment[1] in "AU" else "AU_terminal_on_GC_penultimate"
        )
        correction = corrections[correction_key]
    value = stack_sum + correction
    if not math.isfinite(value):
        raise Stage08Error("non-finite terminal stability")
    return value


def calculate_asymmetry(
    guide_rna: str, passenger_rna: str, terminal_length_nt: int,
    stacks: dict[str, float], corrections: dict[str, float],
) -> tuple[float, float, float]:
    guide = normalize_rna(guide_rna)
    passenger = normalize_rna(passenger_rna)
    if guide != reverse_complement_rna(passenger):
        raise Stage08Error("asymmetry input is not a perfect full-length reverse-complement duplex")
    guide_dg = terminal_dg(guide, terminal_length_nt, stacks, corrections)
    passenger_dg = terminal_dg(passenger, terminal_length_nt, stacks, corrections)
    return guide_dg, passenger_dg, guide_dg - passenger_dg


RNAFOLD_STRUCTURE = re.compile(
    r"^([().]+)\s+\(\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*\)\s*$"
)


def parse_rnafold_output(
    stdout: str, expected_guides: dict[str, str]
) -> dict[str, tuple[float, str]]:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    output: dict[str, tuple[float, str]] = {}
    index = 0
    while index < len(lines):
        if not lines[index].startswith(">") or index + 2 >= len(lines):
            raise Stage08Error("malformed RNAfold multi-FASTA output")
        candidate_id = lines[index][1:].split()[0]
        sequence = normalize_rna(lines[index + 1])
        match = RNAFOLD_STRUCTURE.match(lines[index + 2])
        if not match or candidate_id not in expected_guides:
            raise Stage08Error(f"malformed or unknown RNAfold record: {candidate_id}")
        structure, mfe_text = match.groups()
        mfe = float(mfe_text)
        if sequence != expected_guides[candidate_id]:
            raise Stage08Error(f"RNAfold echoed a different guide sequence for {candidate_id}")
        if len(structure) != len(sequence) or not math.isfinite(mfe):
            raise Stage08Error(f"invalid RNAfold structure/MFE for {candidate_id}")
        output[candidate_id] = (mfe, structure)
        index += 3
    if set(output) != set(expected_guides):
        raise Stage08Error("RNAfold output candidate IDs do not match the candidate table")
    return output


def run_rnafold(
    guides: dict[str, str], temperature_c: float, executable: str = "RNAfold"
) -> dict[str, tuple[float, str]]:
    fasta = "".join(f">{candidate_id}\n{guide}\n" for candidate_id, guide in guides.items())
    completed = subprocess.run(
        [executable, "--noPS", f"--temp={temperature_c}"], input=fasta,
        text=True, capture_output=True, check=False,
    )
    if completed.returncode != 0:
        raise Stage08Error(f"RNAfold failed ({completed.returncode}): {completed.stderr.strip()}")
    return parse_rnafold_output(completed.stdout, guides)


def executable_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable, "--version"], text=True, capture_output=True, check=False
        )
    except OSError as exc:
        raise Stage08Error(f"required executable unavailable: {executable}: {exc}") from exc
    if completed.returncode != 0:
        raise Stage08Error(f"cannot query {executable} version: {completed.stderr.strip()}")
    match = re.search(r"(?:RNAplfold|RNAfold)\s+(\d+\.\d+\.\d+)", completed.stdout + completed.stderr)
    if not match:
        raise Stage08Error(f"cannot parse {executable} version: {(completed.stdout + completed.stderr).strip()}")
    return match.group(1)


def load_stage08_config(path: Path) -> dict[str, object]:
    try:
        configuration = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise Stage08Error(f"cannot parse analysis config {path}: {exc}") from exc
    if "stage08" not in configuration:
        raise Stage08Error("analysis config has no stage08 section")
    stage08 = configuration["stage08"]
    required = {
        "viennarna_version", "temperature_c", "accessibility_parameter_sets",
        "seed_accessibility_guide_positions", "terminal_thermo_model",
        "terminal_thermo_doi", "asymmetry_terminal_lengths_nt",
        "candidate_duplex_overhang_assumption", "self_fold_method",
    }
    if not isinstance(stage08, dict) or required - set(stage08):
        raise Stage08Error("analysis config Stage 08 section is incomplete")
    if stage08["seed_accessibility_guide_positions"] != [2, 8]:
        raise Stage08Error("Stage 08 supports only the canonical g2-g8 seed interval")
    if stage08["asymmetry_terminal_lengths_nt"] != [4, 5]:
        raise Stage08Error("Stage 08 supports only canonical 4-bp primary and 5-bp sensitivity ends")
    if stage08["terminal_thermo_model"] != "Zuber_2022" or stage08["terminal_thermo_doi"] != ZUBER_DOI:
        raise Stage08Error("Stage 08 terminal-thermodynamics configuration is not canonical Zuber 2022")
    if stage08["candidate_duplex_overhang_assumption"] != "none_perfect_complementary_duplex":
        raise Stage08Error("Stage 08 must use a perfect complementary duplex with no overhang")
    return stage08


def _qc(status: str, check: str, value: object, details: str = "", target_id: str = "ALL", length: object = "ALL") -> dict[str, object]:
    return {
        "status": status, "check": check, "target_id": target_id,
        "candidate_length_nt": length, "value": value, "details": details,
    }


def run_stage08(
    candidate_path: Path, manifest_path: Path, zuber_path: Path,
    analysis_config_path: Path, output_root: Path,
) -> dict[str, object]:
    started = time.monotonic()
    stage08 = load_stage08_config(analysis_config_path)
    required_version = str(stage08["viennarna_version"])
    plfold_version = executable_version("RNAplfold")
    rnafold_version = executable_version("RNAfold")
    if plfold_version != required_version or rnafold_version != required_version:
        raise Stage08Error(
            f"ViennaRNA version mismatch: required {required_version}; "
            f"RNAplfold={plfold_version}; RNAfold={rnafold_version}"
        )
    transcripts = load_registered_transcripts(manifest_path)
    candidates = load_and_validate_candidates(candidate_path, transcripts)
    stacks, corrections = load_zuber_parameters(zuber_path)
    accessibility, parameter_rows, fold_runs = calculate_accessibilities(
        candidates, transcripts, list(stage08["accessibility_parameter_sets"]),
        float(stage08["temperature_c"]),
    )
    guides = {
        str(row["candidate_id"]): str(row["antisense_guide_sequence_rna"])
        for row in candidates
    }
    self_folds = run_rnafold(guides, float(stage08["temperature_c"]))

    output_rows: list[dict[str, object]] = []
    for candidate in candidates:
        row = dict(candidate)
        candidate_id = str(candidate["candidate_id"])
        row.update(accessibility[candidate_id])
        guide = str(candidate["antisense_guide_sequence_rna"])
        passenger = str(candidate["target_sequence_rna"])
        for terminal_length in (4, 5):
            guide_dg, passenger_dg, asymmetry = calculate_asymmetry(
                guide, passenger, terminal_length, stacks, corrections
            )
            row[f"guide_5p_terminal_dg_{terminal_length}bp"] = guide_dg
            row[f"passenger_5p_terminal_dg_{terminal_length}bp"] = passenger_dg
            row[f"asymmetry_ddg_{terminal_length}bp"] = asymmetry
        mfe, structure = self_folds[candidate_id]
        row["guide_self_fold_mfe_kcal_mol"] = mfe
        row["guide_self_fold_structure"] = structure
        output_rows.append(row)

    candidate_sha = file_sha256(candidate_path)
    manifest_sha = file_sha256(manifest_path)
    zuber_sha = file_sha256(zuber_path)
    for row in parameter_rows:
        key = (str(row["target_id"]), str(row["transcript_id"]))
        row.update({
            "viennarna_version": required_version,
            "terminal_thermo_parameter_source": f"Zuber_2022_DOI_{ZUBER_DOI}",
            "terminal_thermo_parameter_resource": str(zuber_path),
            "terminal_thermo_parameter_resource_sha256": zuber_sha,
            "candidate_duplex_overhang_assumption": stage08["candidate_duplex_overhang_assumption"],
            "self_fold_method": stage08["self_fold_method"],
            "stage06_candidate_table_sha256": candidate_sha,
            "target_manifest_sha256": manifest_sha,
            "transcript_sequence_sha256_uppercase_dna": transcripts[key]["sequence_sha256"],
        })

    accessibility_fields = [field for pair in ACCESSIBILITY_OUTPUTS.values() for field in pair]
    accessibility_values = [
        float(row[field]) for row in output_rows for field in accessibility_fields
        if row[field] is not None
    ]
    asymmetry_values = [
        float(row[field]) for row in output_rows
        for field in (
            "guide_5p_terminal_dg_4bp", "passenger_5p_terminal_dg_4bp", "asymmetry_ddg_4bp",
            "guide_5p_terminal_dg_5bp", "passenger_5p_terminal_dg_5bp", "asymmetry_ddg_5bp",
        )
    ]
    qc_rows = [
        _qc("PASS", "candidate_count_preserved", len(output_rows), f"input={len(candidates)}"),
        _qc("PASS", "candidate_orientation_and_transcript_slice", len(output_rows)),
        _qc("PASS", "accessibility_probability_range", len(accessibility_values)),
        _qc("PASS", "transcript_level_rnaplfold_reuse", fold_runs,
            f"expected={len(transcripts) * len(stage08['accessibility_parameter_sets'])}"),
        _qc("PASS", "seed_interval_mapping", "[end_1based-7,end_1based-1]"),
        _qc("PASS", "zuber_parameter_integrity", f"{len(stacks)} stacks; {len(corrections)} corrections"),
        _qc("PASS", "terminal_asymmetry_finite", len(asymmetry_values)),
        _qc("PASS", "candidate_duplex_overhang_assumption", stage08["candidate_duplex_overhang_assumption"]),
        _qc("PASS", "guide_self_fold_finite_and_structure_length", len(self_folds)),
        _qc("PASS", "viennarna_version", required_version),
        _qc("PASS", "no_stage07_inputs", 0),
        _qc("PASS", "no_score_rank_weight_gate_or_filter", 0),
    ]
    counts = Counter((str(row["target_id"]), int(row["candidate_length_nt"])) for row in output_rows)
    for (target_id, length), count in sorted(counts.items()):
        qc_rows.append(_qc("INFO", "candidate_count_by_target_length", count, target_id=target_id, length=length))
    if not all(validate_probability(value) == value for value in accessibility_values):
        raise Stage08Error("accessibility range validation failed")
    if not all(math.isfinite(value) for value in asymmetry_values):
        raise Stage08Error("terminal thermodynamics contain non-finite values")
    if fold_runs != len(transcripts) * len(stage08["accessibility_parameter_sets"]):
        raise Stage08Error("RNAplfold transcript-reuse accounting mismatch")

    _write_tsv(output_root / "candidate_biophysics.tsv", BIOPHYSICS_COLUMNS, output_rows)
    _write_tsv(output_root / "stage08_parameters.tsv", PARAMETER_COLUMNS, parameter_rows)
    _write_tsv(output_root / "stage08_qc.tsv", QC_COLUMNS, qc_rows)
    return {
        "runtime_seconds": time.monotonic() - started,
        "rows": len(output_rows), "fold_runs": fold_runs,
        "pass": sum(row["status"] == "PASS" for row in qc_rows),
        "warn": sum(row["status"] == "WARN" for row in qc_rows),
        "fail": sum(row["status"] == "FAIL" for row in qc_rows),
        "info": sum(row["status"] == "INFO" for row in qc_rows),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--target-manifest", required=True, type=Path)
    parser.add_argument("--zuber-resource", required=True, type=Path)
    parser.add_argument("--analysis-config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = run_stage08(
            args.candidates.resolve(), args.target_manifest.resolve(),
            args.zuber_resource.resolve(), args.analysis_config.resolve(),
            args.output_root.resolve(),
        )
    except Stage08Error as exc:
        _write_tsv(
            args.output_root / "stage08_qc.tsv", QC_COLUMNS,
            [_qc("FAIL", "stage08_execution", 1, str(exc))],
        )
        print(f"Stage 08 failed: {exc}")
        return 1
    print(
        f"Stage 08 completed: rows={result['rows']}; "
        f"runtime_seconds={result['runtime_seconds']:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
