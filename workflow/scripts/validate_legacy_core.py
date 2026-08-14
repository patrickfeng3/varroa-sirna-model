#!/usr/bin/env python3
"""Validate the frozen Varroa legacy core without modifying it."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import sys
import time
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


ELIGIBILITY_COLUMNS = [
    "sample", "sample_label", "country", "pooled_mites", "platform",
    "analysis_unit", "biological_virus", "polarity", "exact_mapped_read_names",
    "sense_reads", "antisense_reads", "reference_length_nt", "covered_bases_exact",
    "exact_reference_breadth", "background_usable_bases_depth_masked",
    "background_total_bases", "background_usable_fraction",
    "sample_ambiguous_read_names", "sample_total_competitive_mapped_names",
    "sample_ambiguity_fraction", "pair_ambiguous_read_names_involving_unit",
    "pair_total_assignable_or_ambiguous_names", "pair_ambiguity_fraction",
    "max_strand_total", "largest_length_bin_count", "consensus_last_round",
    "consensus_last_gain_percent", "consensus_last_gain_absolute",
    "consensus_last_accepted_snps", "consensus_converged", "primary_eligible",
    "exploratory_eligible", "failed_primary_checks",
]

READ_FEATURE_COLUMNS = [
    "sample", "mapping_mode", "read_name", "virus", "virus_assignment", "strand",
    "sequence", "length", "five_prime_nt", "three_prime_nt", "count",
    "number_of_mapping_locations", "reference_names",
]

CATALOG_REQUIRED_COLUMNS = {
    "analysis_unit", "biological_virus", "polarity", "seed_id"
}

DEFAULT_SAMPLES = tuple(f"SRR250107{i}" for i in range(50, 71))
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class Check:
    severity: str
    check_id: str
    scope: str
    path: str = ""
    message: str = ""
    observed: str = ""
    expected: str = ""
    size_bytes: str = ""
    sha256: str = ""


@dataclass
class ValidationResult:
    checks: list[Check] = field(default_factory=list)
    required_files: set[Path] = field(default_factory=set)
    started: float = field(default_factory=time.monotonic)
    hashed_bytes: int = 0
    hashed_files: int = 0

    def add(self, severity: str, check_id: str, scope: str, **kwargs) -> None:
        self.checks.append(Check(severity, check_id, scope, **kwargs))

    def require(self, condition: bool, check_id: str, scope: str, **kwargs) -> bool:
        self.add("PASS" if condition else "FAIL", check_id, scope, **kwargs)
        return condition

    @property
    def failed(self) -> bool:
        return any(check.severity == "FAIL" for check in self.checks)


def safe_int(value: str, label: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not an integer: {value!r}") from exc


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def parse_fasta(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    identifier = None
    sequence: list[str] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence).upper()
                identifier = line[1:].split()[0]
                if not identifier:
                    raise ValueError(f"empty FASTA identifier at line {line_number}")
                if identifier in records:
                    raise ValueError(f"duplicate FASTA identifier {identifier!r}")
                sequence = []
            else:
                if identifier is None:
                    raise ValueError(f"sequence before first FASTA header at line {line_number}")
                sequence.append(line)
    if identifier is not None:
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier {identifier!r}")
        records[identifier] = "".join(sequence).upper()
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError("FASTA has no records or contains an empty sequence")
    return records


def reference_span(cigar: str) -> int:
    if cigar == "*":
        raise ValueError("mapped record has '*' CIGAR")
    pieces = CIGAR_RE.findall(cigar)
    if not pieces or "".join(n + op for n, op in pieces) != cigar:
        raise ValueError(f"invalid CIGAR {cigar!r}")
    return sum(int(n) for n, op in pieces if op in "MDN=X")


def path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def check_required_file(result: ValidationResult, root: Path, relative: str) -> Path | None:
    path = root / relative
    exists = path.is_file() and os.access(path, os.R_OK)
    result.require(
        exists, "required_file", relative, path=str(path),
        message="required downstream input is readable" if exists else "required downstream input is missing or unreadable",
    )
    if exists:
        result.required_files.add(path)
        return path
    return None


def validate_feature_table(
    result: ValidationResult,
    path: Path,
    sample: str,
    eligibility_pairs: set[tuple[str, str]],
    expected_exact_assigned_pairs: set[tuple[str, str]],
) -> None:
    categories: dict[str, set[str]] = defaultdict(set)
    rows = 0
    assigned_pairs: set[tuple[str, str]] = set()
    errors: list[str] = []
    try:
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = list(reader.fieldnames or [])
            if header != READ_FEATURE_COLUMNS:
                result.add(
                    "FAIL", "read_features", sample, path=str(path),
                    message=f"ordered schema differs: {header!r}",
                    observed="0", expected=str(READ_FEATURE_COLUMNS),
                )
                result.add(
                    "INFO", "read_feature_categories", sample, path=str(path),
                    message="categories unavailable because the table schema is invalid",
                )
                return
            for line_number, row in enumerate(reader, 2):
                rows += 1
                if None in row:
                    errors.append(f"line {line_number}: surplus fields")
                    break
                if row["sample"] != sample:
                    errors.append(f"line {line_number}: sample {row['sample']!r} != {sample!r}")
                    break
                for column in ("mapping_mode", "virus_assignment", "strand"):
                    categories[column].add(row[column])
                try:
                    length = safe_int(row["length"], "length")
                    count = safe_int(row["count"], "count")
                    loci = safe_int(row["number_of_mapping_locations"], "number_of_mapping_locations")
                    if length != len(row["sequence"]) or count < 0 or loci < 0:
                        raise ValueError("invalid length/count/location relationship")
                except ValueError as exc:
                    errors.append(f"line {line_number}: {exc}")
                    break
                if row["virus_assignment"] == "assigned":
                    pair = (sample, row["virus"])
                    if pair not in eligibility_pairs:
                        errors.append(f"line {line_number}: assigned pair absent from eligibility: {pair}")
                        break
                    if row["mapping_mode"] == "exact":
                        assigned_pairs.add(pair)
    except (OSError, EOFError, UnicodeError, csv.Error) as exc:
        errors.append(f"compressed TSV could not be read completely: {exc}")

    missing_pairs = expected_exact_assigned_pairs - assigned_pairs
    if missing_pairs:
        errors.append(
            "positive exact_mapped_read_names pair(s) lack an exact assigned row: "
            + ",".join(f"{pair[0]}:{pair[1]}" for pair in sorted(missing_pairs))
        )

    result.require(
        not errors, "read_features", sample, path=str(path),
        message="; ".join(errors[:3]) if errors else "gzip stream, schema, rows, and assigned identities are valid",
        observed=str(rows), expected="structurally valid rows",
    )
    result.add(
        "INFO", "read_feature_categories", sample, path=str(path),
        message="observed categories (reported, not constrained by Stage 01 filters)",
        observed="; ".join(f"{key}={sorted(values)}" for key, values in sorted(categories.items())),
    )


def validate_sam(result: ValidationResult, path: Path, sample: str, all_records: OrderedDict[str, str]) -> None:
    sq: OrderedDict[str, int] = OrderedDict()
    pg_lines: list[str] = []
    errors: list[str] = []
    records = 0
    nm_present = 0
    try:
        with path.open(errors="strict") as handle:
            in_records = False
            for line_number, raw in enumerate(handle, 1):
                if raw.startswith("@") and not in_records:
                    line = raw.rstrip("\n")
                    if line.startswith("@SQ\t"):
                        tags = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
                        if "SN" not in tags or "LN" not in tags:
                            errors.append(f"line {line_number}: invalid @SQ")
                            continue
                        try:
                            length = int(tags["LN"])
                        except ValueError:
                            errors.append(f"line {line_number}: non-integer @SQ LN")
                            continue
                        if tags["SN"] in sq or length <= 0:
                            errors.append(f"line {line_number}: duplicate/invalid @SQ")
                        sq[tags["SN"]] = length
                    elif line.startswith("@PG\t"):
                        pg_lines.append(line)
                    continue
                in_records = True
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 11:
                    errors.append(f"line {line_number}: fewer than 11 SAM fields")
                    continue
                records += 1
                try:
                    flag = int(fields[1]); pos = int(fields[3])
                    if flag & 0x4:
                        raise ValueError("unmapped record in --no-unal exact SAM")
                    rname = fields[2]
                    if rname not in sq:
                        raise ValueError(f"RNAME {rname!r} absent from @SQ")
                    span = reference_span(fields[5])
                    if pos < 1 or span < 1 or pos + span - 1 > sq[rname]:
                        raise ValueError(f"alignment outside reference bounds: POS={pos}, span={span}, LN={sq[rname]}")
                    for tag in fields[11:]:
                        if tag.startswith("NM:i:"):
                            nm_present += 1
                            if tag != "NM:i:0":
                                raise ValueError(f"non-zero exact-mapping tag {tag}")
                except ValueError as exc:
                    errors.append(f"line {line_number}: {exc}")
                if len(errors) >= 20:
                    break
    except (OSError, UnicodeError) as exc:
        errors.append(f"SAM could not be read: {exc}")

    fasta_lengths = OrderedDict((name, len(seq)) for name, seq in all_records.items())
    if dict(sq) != dict(fasta_lengths):
        errors.append("@SQ names/lengths differ from current all-virus final FASTA")
    bowtie_exact = any(
        "ID:Bowtie" in line and "VN:1.3.1" in line and re.search(r"(?:^|\s)-v\s+0(?:\s|$)", line)
        for line in pg_lines
    )
    if not bowtie_exact:
        errors.append("@PG does not establish Bowtie 1.3.1 exact mapping with -v 0")
    result.require(
        not errors, "exact_sam", sample, path=str(path),
        message="; ".join(errors[:5]) if errors else "headers, reference dictionary, records, bounds, and exact-mapping provenance are valid",
        observed=f"{records} records; {len(sq)} @SQ; NM present on {nm_present}",
    )
    historic_paths = sorted(set(re.findall(r"/(?:Users|home)/[^\s\"]+", " ".join(pg_lines))))
    result.add(
        "INFO", "sam_pg_historical_paths", sample, path=str(path),
        message="embedded @PG paths are historical diagnostics and are not resolved or used for validation",
        observed=str(len(historic_paths)),
    )


def sha256_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return size, digest.hexdigest()


def validate_core(
    legacy_core: Path,
    output_tsv: Path,
    output_md: Path,
    expected_samples: Iterable[str] = DEFAULT_SAMPLES,
    expected_pair_count: int = 191,
) -> ValidationResult:
    result = ValidationResult()
    root = legacy_core.expanduser().resolve()
    out_tsv = output_tsv.expanduser().resolve()
    out_md = output_md.expanduser().resolve()
    result.require(root.is_dir() and os.access(root, os.R_OK), "legacy_core", "global", path=str(root), message="legacy core is readable")
    output_paths_safe = not path_is_within(out_tsv, root) and not path_is_within(out_md, root)
    result.require(
        output_paths_safe,
        "output_path_safety", "global", observed=f"{out_tsv}; {out_md}", expected="outside legacy core",
        message="all Stage 00 outputs resolve outside the frozen legacy core",
    )

    required_meta = [
        "results/descriptive/eligibility.tsv", "config/virus_catalog.tsv",
        "config/generated_analysis_manifest.tsv", "config/preprocessing_modes.tsv",
        "qc/audit/adapter_audit_summary.tsv",
    ]
    meta = {rel: check_required_file(result, root, rel) for rel in required_meta}
    eligibility_path = meta[required_meta[0]]
    rows: list[dict[str, str]] = []
    if eligibility_path:
        try:
            header, rows = read_tsv(eligibility_path)
            schema_valid = header == ELIGIBILITY_COLUMNS
            result.require(schema_valid, "eligibility_schema", "eligibility", path=str(eligibility_path), observed=str(header), expected=str(ELIGIBILITY_COLUMNS), message="ordered 33-column schema matches")
            if not schema_valid:
                rows = []
        except Exception as exc:
            result.add("FAIL", "eligibility_read", "eligibility", path=str(eligibility_path), message=str(exc))

    expected_sample_set = set(expected_samples)
    samples = sorted({row.get("sample", "") for row in rows if row.get("sample")})
    pairs = {(row.get("sample", ""), row.get("analysis_unit", "")) for row in rows}
    result.require(len(rows) == expected_pair_count, "eligibility_row_count", "eligibility", observed=str(len(rows)), expected=str(expected_pair_count), message="eligibility row count")
    result.require(len(pairs) == len(rows), "eligibility_pair_uniqueness", "eligibility", observed=str(len(pairs)), expected=str(len(rows)), message="sample-analysis_unit pairs are unique")
    result.require(set(samples) == expected_sample_set, "sample_set", "eligibility", observed=",".join(samples), expected=",".join(sorted(expected_sample_set)), message="sample identifiers match the frozen cohort")

    # Required catalogue and manifest consistency.
    if meta["config/virus_catalog.tsv"] and rows:
        try:
            header, catalog = read_tsv(meta["config/virus_catalog.tsv"])
            missing = CATALOG_REQUIRED_COLUMNS - set(header)
            by_unit: dict[str, set[tuple[str, str]]] = defaultdict(set)
            for item in catalog:
                by_unit[item["analysis_unit"]].add((item["biological_virus"], item["polarity"]))
            conflicts = []
            for row in rows:
                expected = (row["biological_virus"], row["polarity"])
                if expected not in by_unit.get(row["analysis_unit"], set()):
                    conflicts.append(f"{row['analysis_unit']}:{expected}")
            result.require(not missing and not conflicts, "virus_catalog_consistency", "catalog", path=str(meta["config/virus_catalog.tsv"]), message=(f"missing columns={sorted(missing)}; conflicts={conflicts[:5]}" if missing or conflicts else "eligibility units, biological viruses, and polarities occur in catalogue"))
        except Exception as exc:
            result.add("FAIL", "virus_catalog_read", "catalog", message=str(exc))

    if meta["config/generated_analysis_manifest.tsv"] and rows:
        try:
            header, manifest = read_tsv(meta["config/generated_analysis_manifest.tsv"])
            manifest_pairs = {(row["sample"], row["virus"]) for row in manifest}
            result.require(header == ["sample", "fastq", "virus", "seed_references"] and manifest_pairs == pairs and len(manifest) == len(pairs), "analysis_manifest_consistency", "manifest", path=str(meta["config/generated_analysis_manifest.tsv"]), observed=f"{len(manifest)} rows", expected=f"{len(pairs)} eligibility pairs", message="selected sample-reference manifest agrees with eligibility")
        except Exception as exc:
            result.add("FAIL", "analysis_manifest_read", "manifest", message=str(exc))

    if meta["config/preprocessing_modes.tsv"]:
        try:
            header, preprocessing = read_tsv(meta["config/preprocessing_modes.tsv"])
            runs = {row.get("run", "") for row in preprocessing}
            result.require({"run", "expected_mode", "evidence"} <= set(header) and runs == expected_sample_set and len(preprocessing) == len(expected_sample_set), "preprocessing_provenance", "preprocessing", path=str(meta["config/preprocessing_modes.tsv"]), message="corrected preprocessing modes/evidence cover the frozen cohort")
        except Exception as exc:
            result.add("FAIL", "preprocessing_provenance_read", "preprocessing", message=str(exc))

    if meta["qc/audit/adapter_audit_summary.tsv"]:
        try:
            header, audits = read_tsv(meta["qc/audit/adapter_audit_summary.tsv"])
            audit_samples = {row.get("sample", "") for row in audits}
            statuses = sorted({row.get("status", "") for row in audits})
            result.require(audit_samples == expected_sample_set and len(audits) == len(expected_sample_set) and statuses == ["PASS"], "preprocessing_audit_summary", "preprocessing", path=str(meta["qc/audit/adapter_audit_summary.tsv"]), observed=str(statuses), expected="['PASS']", message="corrected preprocessing audit covers all samples and passes")
        except Exception as exc:
            result.add("FAIL", "preprocessing_audit_read", "preprocessing", message=str(exc))

    # Concrete downstream inputs.
    feature_paths: dict[str, Path] = {}
    sam_paths: dict[str, Path] = {}
    all_fasta_paths: dict[str, Path] = {}
    pair_fasta_paths: dict[tuple[str, str], Path] = {}
    background_paths: dict[tuple[str, str], Path] = {}
    for sample in samples:
        for relative, target in [
            (f"tables/{sample}/{sample}.read_level_features.tsv.gz", feature_paths),
            (f"alignments/{sample}.all_viruses.exact.sam", sam_paths),
            (f"references/consensus/{sample}.all_viruses.final.fa", all_fasta_paths),
        ]:
            path = check_required_file(result, root, relative)
            if path: target[sample] = path
        audit = root / f"qc/audit/{sample}.adapter_audit.tsv"
        result.add(
            "INFO" if audit.is_file() else "WARN",
            "optional_preprocessing_audit_detail", sample, path=str(audit),
            message="individual historical preprocessing audit is present" if audit.is_file()
            else "individual historical preprocessing audit is absent; required summary/mode provenance remains authoritative",
        )
    for sample, unit in sorted(pairs):
        for suffix, target in [
            ("final.fa", pair_fasta_paths),
            ("final.background_masked.fa", background_paths),
        ]:
            path = check_required_file(result, root, f"references/consensus/{sample}.{unit}.{suffix}")
            if path: target[(sample, unit)] = path

    # FASTA cross-consistency (multi-record safe).
    all_records_by_sample: dict[str, OrderedDict[str, str]] = {}
    for sample, path in all_fasta_paths.items():
        try:
            all_records_by_sample[sample] = parse_fasta(path)
            result.add("PASS", "all_virus_fasta", sample, path=str(path), message="sample-level all-virus FASTA is parseable", observed=str(len(all_records_by_sample[sample])))
        except Exception as exc:
            result.add("FAIL", "all_virus_fasta", sample, path=str(path), message=str(exc))

    rows_by_pair = {(row["sample"], row["analysis_unit"]): row for row in rows}
    pair_ids_by_sample: dict[str, set[str]] = defaultdict(set)
    for pair in sorted(pairs):
        if pair not in pair_fasta_paths or pair not in background_paths or pair[0] not in all_records_by_sample:
            continue
        sample, unit = pair
        try:
            final = parse_fasta(pair_fasta_paths[pair])
            background = parse_fasta(background_paths[pair])
            all_records = all_records_by_sample[sample]
            ids = list(final)
            errors = []
            relevant_all = {
                identifier: sequence
                for identifier, sequence in all_records.items()
                if identifier.split("|", 1)[0] == unit
            }
            wrong_owners = [
                identifier for identifier in ids
                if identifier.split("|", 1)[0] != unit
            ]
            if wrong_owners:
                errors.append(f"record identifiers do not belong to analysis_unit {unit}: {wrong_owners[:5]}")
            if set(background) != set(ids):
                errors.append("background record identifiers differ from pair final FASTA")
            if set(final) != set(relevant_all):
                errors.append("pair final record identifiers differ from the analysis_unit records in all-virus FASTA")
            for identifier in ids:
                seq = final[identifier]
                bg = background.get(identifier, "")
                if relevant_all.get(identifier) != seq:
                    errors.append(f"{identifier}: all-virus sequence differs")
                if len(bg) != len(seq):
                    errors.append(f"{identifier}: background length differs")
                elif any(masked != original and masked != "N" for original, masked in zip(seq, bg)):
                    errors.append(f"{identifier}: background contains a non-N substitution")
            reference_length = sum(map(len, final.values()))
            background_total = sum(map(len, background.values()))
            usable = sum(base != "N" for seq in background.values() for base in seq)
            eligibility = rows_by_pair[pair]
            if reference_length != safe_int(eligibility["reference_length_nt"], "reference_length_nt"):
                errors.append("reference length differs from eligibility")
            if background_total != safe_int(eligibility["background_total_bases"], "background_total_bases"):
                errors.append("background total differs from eligibility")
            if usable != safe_int(eligibility["background_usable_bases_depth_masked"], "background_usable_bases_depth_masked"):
                errors.append("background non-N usable bases differ from eligibility")
            pair_ids_by_sample[sample].update(ids)
            result.require(not errors, "fasta_pair_cross_consistency", f"{sample}:{unit}", path=str(pair_fasta_paths[pair]), message="; ".join(errors[:5]) if errors else "record identifiers, sequences, masking, structure, and eligibility lengths agree", observed=f"{len(ids)} records; reference={reference_length}; background={background_total}; usable={usable}")
        except Exception as exc:
            result.add("FAIL", "fasta_pair_cross_consistency", f"{sample}:{unit}", message=str(exc))
    for sample, all_records in all_records_by_sample.items():
        result.require(set(all_records) == pair_ids_by_sample[sample], "all_virus_fasta_record_set", sample, path=str(all_fasta_paths[sample]), observed=str(len(all_records)), expected=str(len(pair_ids_by_sample[sample])), message="all-virus FASTA is exactly the union of selected pair FASTA records")

    positive_exact_pairs = {
        (row["sample"], row["analysis_unit"])
        for row in rows
        if safe_int(row["exact_mapped_read_names"], "exact_mapped_read_names") > 0
    }
    for sample, path in feature_paths.items():
        validate_feature_table(
            result, path, sample, pairs,
            {pair for pair in positive_exact_pairs if pair[0] == sample},
        )
    for sample, path in sam_paths.items():
        if sample in all_records_by_sample:
            validate_sam(result, path, sample, all_records_by_sample[sample])

    # Optional/historical diagnostics never fail Stage 00.
    for relative in ["tables/all_samples.mapping_summary.tsv", "results/provenance/file_manifest.tsv", "MANIFEST.sha256"]:
        path = root / relative
        result.add("INFO" if path.is_file() else "WARN", "optional_historical_artifact", relative, path=str(path), message="historical diagnostic is present" if path.is_file() else "historical diagnostic is absent; canonical downstream interpretation is unaffected")
    legacy_manifest = root / "results/provenance/file_manifest.tsv"
    if legacy_manifest.is_file():
        try:
            _, manifest_rows = read_tsv(legacy_manifest)
            covered = {root / row["path"] for row in manifest_rows}
            absent = len(result.required_files - covered)
            result.add("WARN" if absent else "INFO", "legacy_identity_coverage", "provenance", path=str(legacy_manifest), observed=str(absent), expected="0", message="required files absent from the historical identity manifest; Stage 00 records fresh identities" if absent else "historical identity manifest covers required files")
        except Exception as exc:
            result.add("WARN", "legacy_identity_manifest_read", "provenance", path=str(legacy_manifest), message=str(exc))

    # Hash every unique required file exactly once.
    for path in sorted(result.required_files):
        relative = str(path.relative_to(root))
        try:
            size, digest = sha256_file(path)
            result.hashed_bytes += size
            result.hashed_files += 1
            result.add("INFO", "required_file_identity", relative, path=str(path), message="fresh frozen-input identity", size_bytes=str(size), sha256=digest)
        except OSError as exc:
            result.add("FAIL", "required_file_identity", relative, path=str(path), message=f"hashing failed: {exc}")

    elapsed = time.monotonic() - result.started
    result.add("INFO", "validation_runtime", "global", observed=f"{elapsed:.3f}", expected="seconds", message="wall-clock runtime measured inside validator")
    result.add("INFO", "hash_volume", "global", observed=str(result.hashed_bytes), expected="bytes", message=f"{result.hashed_files} unique required files hashed once")
    if output_paths_safe:
        write_reports(result, out_tsv, out_md)
    return result


def write_reports(result: ValidationResult, output_tsv: Path, output_md: Path) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    fields = list(Check.__dataclass_fields__)
    tsv_tmp = output_tsv.with_suffix(output_tsv.suffix + ".tmp")
    with tsv_tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for check in result.checks:
            writer.writerow(vars(check))
    os.replace(tsv_tmp, output_tsv)

    counts = Counter(check.severity for check in result.checks)
    overall = "FAIL" if result.failed else "PASS"
    noteworthy = [check for check in result.checks if check.severity in {"FAIL", "WARN"}]
    md_tmp = output_md.with_suffix(output_md.suffix + ".tmp")
    with md_tmp.open("w") as handle:
        handle.write("# Legacy core validation\n\n")
        handle.write(f"**Overall status:** {overall}\n\n")
        handle.write(f"Checks: {counts['PASS']} PASS, {counts['FAIL']} FAIL, {counts['WARN']} WARN, {counts['INFO']} INFO.\n\n")
        handle.write(f"Required-file identities: {result.hashed_files} files, {result.hashed_bytes} bytes SHA-256 hashed.\n\n")
        handle.write("## FAIL and WARN checks\n\n")
        if not noteworthy:
            handle.write("None.\n")
        else:
            for check in noteworthy:
                handle.write(f"- **{check.severity} `{check.check_id}`** ({check.scope}): {check.message}\n")
        handle.write("\nFull machine-readable details, including per-file identities, are in `legacy_core_validation.tsv`.\n")
    os.replace(md_tmp, output_md)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-core", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-md", required=True, type=Path)
    args = parser.parse_args(argv)
    result = validate_core(args.legacy_core, args.output_tsv, args.output_md)
    return 1 if result.failed else 0


if __name__ == "__main__":
    sys.exit(main())
