STAGE01_ROOT = "results/01_viral_23_24"


rule stage01_viral_23_24:
    input:
        validator="results/00_validation/legacy_core_validation.tsv",
        script="workflow/scripts/stage01.py",
        analysis_config="config/analysis.yaml",
    output:
        accounting=f"{STAGE01_ROOT}/qc/stage01_accounting.tsv",
        length_pair=f"{STAGE01_ROOT}/length_spectrum/length_distribution_by_pair.tsv",
        length_sample=f"{STAGE01_ROOT}/length_spectrum/length_distribution_by_sample.tsv",
        length_dataset=f"{STAGE01_ROOT}/length_spectrum/length_distribution_across_dataset.tsv",
        counts_pair=f"{STAGE01_ROOT}/fixed_23_24/23_24_counts_by_pair.tsv",
        fractions_pair=f"{STAGE01_ROOT}/fixed_23_24/23_24_fractions_by_pair.tsv",
        fixed_sample=f"{STAGE01_ROOT}/fixed_23_24/23_24_by_sample.tsv",
        fixed_dataset=f"{STAGE01_ROOT}/fixed_23_24/23_24_across_dataset.tsv",
    params:
        legacy_core=str(LEGACY_CORE),
        output_root=STAGE01_ROOT,
    shell:
        "python3 {input.script:q} --legacy-core {params.legacy_core:q} "
        "--config {input.analysis_config:q} --output-root {params.output_root:q}"
