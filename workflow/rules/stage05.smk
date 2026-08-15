STAGE05_ROOT = "results/05_viral_transitivity"


rule stage05_viral_transitivity:
    input:
        script="workflow/scripts/stage05.py",
        analysis_config="config/stage05.yaml",
    output:
        coordinate_qc=f"{STAGE05_ROOT}/coordinate_qc.tsv",
        eligible=f"{STAGE05_ROOT}/eligible_positive_sense_units.tsv",
        historical_pair=f"{STAGE05_ROOT}/historical_v1.4.1_replication/transitivity_by_pair.tsv",
        historical_pair_results=f"{STAGE05_ROOT}/historical_v1.4.1_replication/pair_balanced_results.tsv",
        historical_virus_results=f"{STAGE05_ROOT}/historical_v1.4.1_replication/virus_balanced_results.tsv",
        historical_loo=f"{STAGE05_ROOT}/historical_v1.4.1_replication/leave_one_virus_out.tsv",
        crosscorr=f"{STAGE05_ROOT}/historical_v1.4.1_replication/cross_correlation.tsv",
        regression=f"{STAGE05_ROOT}/historical_v1.4.1_replication/regression_check.tsv",
        canonical_pair=f"{STAGE05_ROOT}/canonical_transitivity_analysis/transitivity_by_pair.tsv",
        canonical_sample=f"{STAGE05_ROOT}/canonical_transitivity_analysis/transitivity_by_sample.tsv",
        canonical_results=f"{STAGE05_ROOT}/canonical_transitivity_analysis/sample_balanced_results.tsv",
        pair_sensitivity=f"{STAGE05_ROOT}/canonical_transitivity_analysis/pair_balanced_sensitivity.tsv",
        virus_sensitivity=f"{STAGE05_ROOT}/canonical_transitivity_analysis/virus_balanced_sensitivity.tsv",
        canonical_loo=f"{STAGE05_ROOT}/canonical_transitivity_analysis/leave_one_virus_out.tsv",
        multiple_testing=f"{STAGE05_ROOT}/canonical_transitivity_analysis/multiple_testing_summary.tsv",
        final=f"{STAGE05_ROOT}/canonical_transitivity_analysis/final_transitivity_summary.tsv",
    params:
        legacy_core=str(LEGACY_CORE),
        output_root=STAGE05_ROOT,
    shell:
        "python3 {input.script:q} --legacy-core {params.legacy_core:q} "
        "--config {input.analysis_config:q} --output-root {params.output_root:q}"
