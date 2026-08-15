from pathlib import Path


local_paths = Path("config/paths.local.yaml")
if local_paths.exists():
    configfile: str(local_paths)
else:
    configfile: "config/paths.example.yaml"


include: "workflow/rules/validate_legacy_core.smk"
include: "workflow/rules/stage01.smk"
include: "workflow/rules/stage02.smk"
include: "workflow/rules/stage03.smk"
include: "workflow/rules/stage04.smk"
include: "workflow/rules/stage05.smk"


rule all:
    input:
        "results/00_validation/legacy_core_validation.tsv",
        "results/00_validation/legacy_core_validation.md",
        "results/01_viral_23_24/qc/stage01_accounting.tsv",
        "results/01_viral_23_24/length_spectrum/length_distribution_by_pair.tsv",
        "results/01_viral_23_24/length_spectrum/length_distribution_by_sample.tsv",
        "results/01_viral_23_24/length_spectrum/length_distribution_across_dataset.tsv",
        "results/01_viral_23_24/fixed_23_24/23_24_counts_by_pair.tsv",
        "results/01_viral_23_24/fixed_23_24/23_24_fractions_by_pair.tsv",
        "results/01_viral_23_24/fixed_23_24/23_24_by_sample.tsv",
        "results/01_viral_23_24/fixed_23_24/23_24_across_dataset.tsv",
        "results/02_terminal_enrichment/qc/stage02_accounting.tsv",
        "results/02_terminal_enrichment/observed/terminal_observed_by_pair.tsv",
        "results/02_terminal_enrichment/background/terminal_expected_by_pair.tsv",
        "results/02_terminal_enrichment/enrichment/terminal_enrichment_by_pair.tsv",
        "results/02_terminal_enrichment/enrichment/terminal_enrichment_by_sample.tsv",
        "results/02_terminal_enrichment/enrichment/terminal_enrichment_across_dataset.tsv",
        "results/02_terminal_enrichment/enrichment/terminal_enrichment_pooled_abundance.tsv",
        "results/02_terminal_enrichment/comparisons/enrichment_23_vs_24.tsv",
        "results/03_steprna/qc/stage03_accounting.tsv",
        "results/03_steprna/provenance/software_versions.tsv",
        "results/03_steprna/provenance/run_manifest.tsv",
        "results/03_steprna/inputs/input_manifest.tsv",
        "results/03_steprna/inputs/focal_reference_manifest.tsv.gz",
        "results/03_steprna/inputs/passenger_manifest.tsv.gz",
        "results/03_steprna/parsed/passenger_recovery_by_pair.tsv",
        "results/03_steprna/parsed/overhang_spectrum_by_pair.tsv",
        "results/03_steprna/parsed/passenger_length_by_pair.tsv",
        "results/03_steprna/parsed/joint_geometry_by_pair.tsv",
        "results/03_steprna/parsed/joint_geometry_references.tsv.gz",
        "results/03_steprna/parsed/joint_geometry_spectrum_by_pair.tsv",
        "results/03_steprna/qc/stage03_joint_geometry_spectrum_accounting.tsv",
        "results/04_duplex_geometry/qc/stage04_accounting.tsv",
        "results/04_duplex_geometry/population/full_spectrum_by_sample.tsv",
        "results/04_duplex_geometry/population/full_spectrum_across_dataset.tsv",
        "results/04_duplex_geometry/population/passenger_recovery_across_dataset.tsv",
        "results/04_duplex_geometry/population/joint_geometry_by_sample.tsv",
        "results/04_duplex_geometry/population/joint_geometry_across_dataset.tsv",
        "results/04_duplex_geometry/comparisons/paired_23_vs_24.tsv",
        "results/04_duplex_geometry/sequence_features/geometry_terminal_by_pair.tsv",
        "results/04_duplex_geometry/sequence_features/geometry_terminal_by_sample.tsv",
        "results/04_duplex_geometry/sequence_features/geometry_terminal_across_dataset.tsv",
        "results/04_duplex_geometry/sequence_features/geometry_specific_contrasts.tsv",
        "results/04_duplex_geometry/sequence_features/redundancy.tsv",
        "results/04_duplex_geometry/population/joint_geometry_spectrum_by_sample.tsv",
        "results/04_duplex_geometry/population/joint_geometry_spectrum_across_dataset.tsv",
        "results/04_duplex_geometry/population/joint_geometry_mode_by_pair.tsv",
        "results/04_duplex_geometry/population/joint_geometry_spectrum_summary.tsv",
        "results/04_duplex_geometry/qc/stage04_joint_geometry_spectrum_accounting.tsv",
        "results/05_viral_transitivity/coordinate_qc.tsv",
        "results/05_viral_transitivity/eligible_positive_sense_units.tsv",
        "results/05_viral_transitivity/historical_v1.4.1_replication/transitivity_by_pair.tsv",
        "results/05_viral_transitivity/historical_v1.4.1_replication/pair_balanced_results.tsv",
        "results/05_viral_transitivity/historical_v1.4.1_replication/virus_balanced_results.tsv",
        "results/05_viral_transitivity/historical_v1.4.1_replication/leave_one_virus_out.tsv",
        "results/05_viral_transitivity/historical_v1.4.1_replication/cross_correlation.tsv",
        "results/05_viral_transitivity/historical_v1.4.1_replication/regression_check.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/transitivity_by_pair.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/transitivity_by_sample.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/sample_balanced_results.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/pair_balanced_sensitivity.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/virus_balanced_sensitivity.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/leave_one_virus_out.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/multiple_testing_summary.tsv",
        "results/05_viral_transitivity/canonical_transitivity_analysis/final_transitivity_summary.tsv",
