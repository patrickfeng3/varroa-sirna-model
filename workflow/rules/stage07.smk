STAGE07_ROOT = "results/07_empirical_sequence"
STAGE02_REGRESSION_ROOT = "results/02_terminal_enrichment"


rule stage07_empirical_sequence:
    input:
        script="workflow/scripts/stage07.py",
        eligibility=str(LEGACY_CORE / "results/descriptive/eligibility.tsv"),
    output:
        positional_pair=f"{STAGE07_ROOT}/positional_by_pair.tsv",
        positional_sample=f"{STAGE07_ROOT}/positional_by_sample.tsv",
        positional_summary=f"{STAGE07_ROOT}/positional_summary.tsv",
        gc_pair=f"{STAGE07_ROOT}/gc9_14_by_pair.tsv",
        gc_sample=f"{STAGE07_ROOT}/gc9_14_by_sample.tsv",
        gc_summary=f"{STAGE07_ROOT}/gc9_14_summary.tsv",
        regional_pair=f"{STAGE07_ROOT}/regional_gc6_by_pair.tsv",
        regional_sample=f"{STAGE07_ROOT}/regional_gc6_by_sample.tsv",
        regional_summary=f"{STAGE07_ROOT}/regional_gc6_summary.tsv",
        regional_discovery=f"{STAGE07_ROOT}/regional_gc6_discovery.tsv",
        literature=f"{STAGE07_ROOT}/literature_validation.tsv",
        discovery=f"{STAGE07_ROOT}/discovery_summary.tsv",
        sense=f"{STAGE07_ROOT}/sense_comparator.tsv",
        accounting=f"{STAGE07_ROOT}/qc/stage07_accounting.tsv",
        regression=f"{STAGE07_ROOT}/qc/stage02_terminal_regression.tsv",
        provenance=f"{STAGE07_ROOT}/provenance/stage07_manifest.tsv",
    params:
        legacy_core=str(LEGACY_CORE),
        stage02_pair_reference=f"{STAGE02_REGRESSION_ROOT}/enrichment/terminal_enrichment_by_pair.tsv",
        stage02_across_reference=f"{STAGE02_REGRESSION_ROOT}/enrichment/terminal_enrichment_across_dataset.tsv",
        output_root=STAGE07_ROOT,
    shell:
        "python3 {input.script:q} --legacy-core {params.legacy_core:q} "
        "--stage02-pair-reference {params.stage02_pair_reference:q} "
        "--stage02-across-reference {params.stage02_across_reference:q} "
        "--output-root {params.output_root:q}"
