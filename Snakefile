from pathlib import Path


local_paths = Path("config/paths.local.yaml")
if local_paths.exists():
    configfile: str(local_paths)
else:
    configfile: "config/paths.example.yaml"


include: "workflow/rules/validate_legacy_core.smk"
include: "workflow/rules/stage01.smk"
include: "workflow/rules/stage02.smk"


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
