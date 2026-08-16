STAGE10_ROOT = "results/10_candidate_integration"


rule stage10a_candidate_integration:
    input:
        script="workflow/scripts/stage10a.py",
    output:
        candidates=f"{STAGE10_ROOT}/candidate_stage10a.tsv",
        correlations=f"{STAGE10_ROOT}/stage10a_layer_correlations.tsv",
        pareto_summary=f"{STAGE10_ROOT}/stage10a_pareto_summary.tsv",
        parameters=f"{STAGE10_ROOT}/stage10_parameters.tsv",
        qc=f"{STAGE10_ROOT}/stage10_qc.tsv",
    params:
        layer1="results/09_feature_layers/09A_layer1_accumulation/candidate_layer1.tsv",
        layer2="results/09_feature_layers/09B_layer2_guide_competence/candidate_layer2.tsv",
        layer3="results/09_feature_layers/09C_layer3_target_engagement/candidate_layer3.tsv",
        output_root=STAGE10_ROOT,
    shell:
        "python3 {input.script:q} --layer1 {params.layer1:q} --layer2 {params.layer2:q} "
        "--layer3 {params.layer3:q} --output-root {params.output_root:q}"
