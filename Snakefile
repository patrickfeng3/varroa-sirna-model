from pathlib import Path


local_paths = Path("config/paths.local.yaml")
if local_paths.exists():
    configfile: str(local_paths)
else:
    configfile: "config/paths.example.yaml"


include: "workflow/rules/validate_legacy_core.smk"
include: "workflow/rules/stage01.smk"


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
