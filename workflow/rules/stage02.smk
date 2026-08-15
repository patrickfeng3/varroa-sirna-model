STAGE02_ROOT = "results/02_terminal_enrichment"


rule stage02_terminal_enrichment:
    input:
        validator="results/00_validation/legacy_core_validation.tsv",
        script="workflow/scripts/stage02.py",
        analysis_config="config/analysis.yaml",
    output:
        accounting=f"{STAGE02_ROOT}/qc/stage02_accounting.tsv",
        observed=f"{STAGE02_ROOT}/observed/terminal_observed_by_pair.tsv",
        expected=f"{STAGE02_ROOT}/background/terminal_expected_by_pair.tsv",
        enrichment_pair=f"{STAGE02_ROOT}/enrichment/terminal_enrichment_by_pair.tsv",
        enrichment_sample=f"{STAGE02_ROOT}/enrichment/terminal_enrichment_by_sample.tsv",
        enrichment_dataset=f"{STAGE02_ROOT}/enrichment/terminal_enrichment_across_dataset.tsv",
        pooled=f"{STAGE02_ROOT}/enrichment/terminal_enrichment_pooled_abundance.tsv",
        comparison=f"{STAGE02_ROOT}/comparisons/enrichment_23_vs_24.tsv",
    params:
        legacy_core=str(LEGACY_CORE),
        output_root=STAGE02_ROOT,
    shell:
        "python3 {input.script:q} --legacy-core {params.legacy_core:q} "
        "--config {input.analysis_config:q} --output-root {params.output_root:q}"
