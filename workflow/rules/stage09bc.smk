STAGE09_ROOT = "results/09_feature_layers"
STAGE09B_ROOT = f"{STAGE09_ROOT}/09B_layer2_guide_competence"
STAGE09C_ROOT = f"{STAGE09_ROOT}/09C_layer3_target_engagement"


rule stage09bc_candidate_layers:
    input:
        script="workflow/scripts/stage09bc.py",
    output:
        layer2=f"{STAGE09B_ROOT}/candidate_layer2.tsv",
        layer2_sensitivity_23=f"{STAGE09B_ROOT}/layer2_weight_sensitivity_23nt.tsv",
        layer2_sensitivity_24=f"{STAGE09B_ROOT}/layer2_weight_sensitivity_24nt.tsv",
        layer2_correlations=f"{STAGE09B_ROOT}/layer2_correlations.tsv",
        layer3=f"{STAGE09C_ROOT}/candidate_layer3.tsv",
        layer3_sensitivity_23=f"{STAGE09C_ROOT}/layer3_weight_sensitivity_23nt.tsv",
        layer3_sensitivity_24=f"{STAGE09C_ROOT}/layer3_weight_sensitivity_24nt.tsv",
        layer3_correlations=f"{STAGE09C_ROOT}/layer3_correlations.tsv",
    params:
        stage08_candidates="results/08_candidate_biophysics/candidate_biophysics.tsv",
        output_root=STAGE09_ROOT,
    shell:
        "python3 {input.script:q} --stage08-candidates {params.stage08_candidates:q} "
        "--output-root {params.output_root:q}"
