import csv
import gzip
import hashlib
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow/scripts/validate_legacy_core.py"
SPEC = importlib.util.spec_from_file_location("validate_legacy_core", SCRIPT)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
sys.modules[SPEC.name] = validator
SPEC.loader.exec_module(validator)


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


class Stage00Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.core = self.base / "legacy"
        self.out = self.base / "repo/results/00_validation"
        self.sample = "SRR1"
        self.unit = "VIRUS"
        self.identifier = "VIRUS|segment1|ref"
        self.sequence = "ACGTACGT"
        self._make_valid_fixture()

    def tearDown(self):
        self.temp.cleanup()

    def _eligibility_row(self):
        values = {column: "0" for column in validator.ELIGIBILITY_COLUMNS}
        values.update({
            "sample": self.sample, "sample_label": "S1", "country": "X",
            "platform": "ILLUMINA", "analysis_unit": self.unit,
            "biological_virus": self.unit, "polarity": "+ssRNA",
            "exact_mapped_read_names": "1",
            "reference_length_nt": "8", "background_total_bases": "8",
            "background_usable_bases_depth_masked": "6",
            "primary_eligible": "true", "exploratory_eligible": "true",
        })
        return [values[column] for column in validator.ELIGIBILITY_COLUMNS]

    def _make_valid_fixture(self):
        c = self.core
        write_tsv(c / "results/descriptive/eligibility.tsv", validator.ELIGIBILITY_COLUMNS, [self._eligibility_row()])
        write_tsv(c / "config/virus_catalog.tsv", ["analysis_unit", "biological_virus", "seed_id", "polarity"], [[self.unit, self.unit, "seed", "+ssRNA"]])
        write_tsv(c / "config/generated_analysis_manifest.tsv", ["sample", "fastq", "virus", "seed_references"], [[self.sample, "raw.fastq.gz", self.unit, "seed.fa"]])
        write_tsv(c / "config/preprocessing_modes.tsv", ["run", "expected_mode", "evidence"], [[self.sample, "already_trimmed", "validated"]])
        write_tsv(c / "qc/audit/adapter_audit_summary.tsv", ["sample", "status"], [[self.sample, "PASS"]])
        write_tsv(c / f"qc/audit/{self.sample}.adapter_audit.tsv", ["sample", "status"], [[self.sample, "PASS"]])
        feature = c / f"tables/{self.sample}/{self.sample}.read_level_features.tsv.gz"
        feature.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(feature, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(validator.READ_FEATURE_COLUMNS)
            writer.writerow([self.sample, "exact", "r1", self.unit, "assigned", "sense", "ACGT", 4, "A", "T", 1, 1, self.identifier])
        consensus = c / "references/consensus"
        consensus.mkdir(parents=True, exist_ok=True)
        (consensus / f"{self.sample}.{self.unit}.final.fa").write_text(f">{self.identifier}\n{self.sequence}\n")
        (consensus / f"{self.sample}.{self.unit}.final.background_masked.fa").write_text(f">{self.identifier}\n{self.sequence[:2]}NN{self.sequence[4:]}\n")
        (consensus / f"{self.sample}.all_viruses.final.fa").write_text(f">{self.identifier}\n{self.sequence}\n")
        sam = c / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.parent.mkdir(parents=True, exist_ok=True)
        sam.write_text(
            f"@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:{self.identifier}\tLN:8\n"
            f"@PG\tID:Bowtie\tVN:1.3.1\tCL:bowtie -v 0 -a --best --strata --sam --no-unal /old/path\n"
            f"r1\t0\t{self.identifier}\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
        )

    def run_validation(self):
        return validator.validate_core(
            self.core, self.out / "legacy_core_validation.tsv",
            self.out / "legacy_core_validation.md",
            expected_samples=[self.sample], expected_pair_count=1,
        )

    def write_feature(self, header, row):
        feature = self.core / f"tables/{self.sample}/{self.sample}.read_level_features.tsv.gz"
        with gzip.open(feature, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerow(row)

    def test_valid_fixture_passes_with_warn_and_info_diagnostics(self):
        result = self.run_validation()
        self.assertFalse(result.failed)
        severities = {check.severity for check in result.checks}
        self.assertIn("PASS", severities)
        self.assertIn("WARN", severities)
        self.assertIn("INFO", severities)
        self.assertEqual(result.hashed_files, 10)

    def test_required_missing_file_fails(self):
        (self.core / f"alignments/{self.sample}.all_viruses.exact.sam").unlink()
        result = self.run_validation()
        self.assertTrue(result.failed)
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "required_file" for c in result.checks))

    def test_missing_feature_column_is_structured_fail_with_reports(self):
        header = validator.READ_FEATURE_COLUMNS[:-1]
        row = [self.sample, "exact", "r1", self.unit, "assigned", "sense", "ACGT", 4, "A", "T", 1, 1]
        self.write_feature(header, row)
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "read_features" for c in result.checks))
        self.assertTrue((self.out / "legacy_core_validation.tsv").is_file())
        self.assertTrue((self.out / "legacy_core_validation.md").is_file())

    def test_reordered_feature_schema_is_structured_fail_with_reports(self):
        header = list(validator.READ_FEATURE_COLUMNS)
        header[0], header[1] = header[1], header[0]
        row = ["exact", self.sample, "r1", self.unit, "assigned", "sense", "ACGT", 4, "A", "T", 1, 1, self.identifier]
        self.write_feature(header, row)
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "read_features" for c in result.checks))
        self.assertTrue((self.out / "legacy_core_validation.tsv").is_file())

    def test_corrupt_gzip_is_structured_fail_with_reports(self):
        feature = self.core / f"tables/{self.sample}/{self.sample}.read_level_features.tsv.gz"
        feature.write_bytes(b"\x1f\x8b\x08corrupt")
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "read_features" for c in result.checks))
        self.assertTrue((self.out / "legacy_core_validation.tsv").is_file())

    def test_positive_eligibility_pair_requires_exact_assigned_row(self):
        row = [self.sample, "exact", "r1", "AMBIGUOUS_MULTI_VIRUS", "ambiguous_multi_virus", "sense", "ACGT", 4, "A", "T", 1, 1, self.identifier]
        self.write_feature(validator.READ_FEATURE_COLUMNS, row)
        result = self.run_validation()
        self.assertTrue(any(
            c.severity == "FAIL" and "lack an exact assigned row" in c.message
            for c in result.checks
        ))

    def test_optional_historical_artifact_is_warn_only(self):
        result = self.run_validation()
        optional = [c for c in result.checks if c.check_id == "optional_historical_artifact"]
        self.assertTrue(optional)
        self.assertTrue(all(c.severity == "WARN" for c in optional))
        self.assertFalse(result.failed)

    def test_missing_individual_adapter_audit_is_warn_only(self):
        (self.core / f"qc/audit/{self.sample}.adapter_audit.tsv").unlink()
        result = self.run_validation()
        checks = [c for c in result.checks if c.check_id == "optional_preprocessing_audit_detail"]
        self.assertEqual([c.severity for c in checks], ["WARN"])
        self.assertFalse(result.failed)

    def test_fasta_all_virus_sequence_mismatch_fails(self):
        path = self.core / f"references/consensus/{self.sample}.all_viruses.final.fa"
        path.write_text(f">{self.identifier}\nACGTTCGT\n")
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "fasta_pair_cross_consistency" for c in result.checks))

    def test_background_non_n_substitution_fails(self):
        path = self.core / f"references/consensus/{self.sample}.{self.unit}.final.background_masked.fa"
        path.write_text(f">{self.identifier}\nACTTNNGT\n")
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "non-N substitution" in c.message for c in result.checks))

    def test_wrong_fasta_analysis_unit_ownership_fails(self):
        wrong = "OTHER|segment1|ref"
        consensus = self.core / "references/consensus"
        (consensus / f"{self.sample}.{self.unit}.final.fa").write_text(f">{wrong}\n{self.sequence}\n")
        (consensus / f"{self.sample}.{self.unit}.final.background_masked.fa").write_text(f">{wrong}\nACNNACGT\n")
        (consensus / f"{self.sample}.all_viruses.final.fa").write_text(f">{wrong}\n{self.sequence}\n")
        sam = self.core / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.write_text(
            f"@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:{wrong}\tLN:8\n"
            f"@PG\tID:Bowtie\tVN:1.3.1\tCL:bowtie -v 0 -a --best --strata --sam --no-unal /old/path\n"
            f"r1\t0\t{wrong}\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
        )
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "do not belong" in c.message for c in result.checks))

    def test_background_length_mismatch_fails(self):
        path = self.core / f"references/consensus/{self.sample}.{self.unit}.final.background_masked.fa"
        path.write_text(f">{self.identifier}\nACNNACG\n")
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "background length differs" in c.message for c in result.checks))

    def test_background_usable_base_mismatch_fails(self):
        path = self.core / f"references/consensus/{self.sample}.{self.unit}.final.background_masked.fa"
        path.write_text(f">{self.identifier}\nNNNNACGT\n")
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "usable bases differ" in c.message for c in result.checks))

    def test_multirecord_fasta_with_reordered_sam_sq_passes(self):
        second = "VIRUS|segment2|ref"
        pair = self.core / f"references/consensus/{self.sample}.{self.unit}.final.fa"
        bg = self.core / f"references/consensus/{self.sample}.{self.unit}.final.background_masked.fa"
        all_fa = self.core / f"references/consensus/{self.sample}.all_viruses.final.fa"
        pair.write_text(f">{self.identifier}\nACGT\n>{second}\nACGT\n")
        bg.write_text(f">{self.identifier}\nACNN\n>{second}\nACGT\n")
        all_fa.write_text(f">{self.identifier}\nACGT\n>{second}\nACGT\n")
        sam = self.core / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.write_text(
            f"@HD\tVN:1.0\tSO:unsorted\n"
            f"@SQ\tSN:{second}\tLN:4\n@SQ\tSN:{self.identifier}\tLN:4\n"
            f"@PG\tID:Bowtie\tVN:1.3.1\tCL:bowtie -v 0 -a --best --strata --sam --no-unal /old/path\n"
            f"r1\t0\t{self.identifier}\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
        )
        result = self.run_validation()
        self.assertFalse(result.failed)

    def test_unknown_sam_rname_fails(self):
        sam = self.core / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.write_text(sam.read_text().replace(f"r1\t0\t{self.identifier}\t", "r1\t0\tUNKNOWN\t"))
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "exact_sam" for c in result.checks))

    def test_incorrect_sam_sq_length_fails(self):
        sam = self.core / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.write_text(sam.read_text().replace("\tLN:8", "\tLN:9"))
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "@SQ names/lengths differ" in c.message for c in result.checks))

    def test_invalid_sam_pos_fails(self):
        sam = self.core / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.write_text(sam.read_text().replace(f"{self.identifier}\t1\t255", f"{self.identifier}\t0\t255"))
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "outside reference bounds" in c.message for c in result.checks))

    def test_nonzero_nm_fails(self):
        sam = self.core / f"alignments/{self.sample}.all_viruses.exact.sam"
        sam.write_text(sam.read_text().replace("NM:i:0", "NM:i:1"))
        result = self.run_validation()
        self.assertTrue(any(c.severity == "FAIL" and "non-zero exact-mapping tag" in c.message for c in result.checks))

    def test_known_sha256_and_size_are_recorded(self):
        target = self.core / "config/virus_catalog.tsv"
        expected_bytes = target.read_bytes()
        result = self.run_validation()
        identity = next(
            c for c in result.checks
            if c.check_id == "required_file_identity" and c.scope == "config/virus_catalog.tsv"
        )
        self.assertEqual(identity.size_bytes, str(len(expected_bytes)))
        self.assertEqual(identity.sha256, hashlib.sha256(expected_bytes).hexdigest())

    def test_cli_failure_returns_nonzero_and_writes_reports(self):
        (self.core / f"alignments/{self.sample}.all_viruses.exact.sam").unlink()
        tsv = self.out / "cli.tsv"
        md = self.out / "cli.md"
        exit_code = validator.main([
            "--legacy-core", str(self.core), "--output-tsv", str(tsv), "--output-md", str(md)
        ])
        self.assertEqual(exit_code, 1)
        self.assertTrue(tsv.is_file())
        self.assertTrue(md.is_file())
        self.assertIn("**Overall status:** FAIL", md.read_text())

    def test_output_inside_legacy_core_fails_safety_check(self):
        result = validator.validate_core(
            self.core, self.core / "results/00_validation/x.tsv",
            self.core / "results/00_validation/x.md",
            expected_samples=[self.sample], expected_pair_count=1,
        )
        self.assertTrue(any(c.severity == "FAIL" and c.check_id == "output_path_safety" for c in result.checks))
        self.assertFalse((self.core / "results/00_validation/x.tsv").exists())
        self.assertFalse((self.core / "results/00_validation/x.md").exists())


if __name__ == "__main__":
    unittest.main()
