STAGE09A_ROOT = "results/09_feature_layers/09A_layer1_accumulation"


rule stage09a_layer1_accumulation:
    input:
        script="workflow/scripts/stage09a.py",
        model_script="workflow/scripts/stage09a_model.R",
        solver="workflow/scripts/stage09a_solver.R",
    output:
        coefficients=f"{STAGE09A_ROOT}/layer1_model_coefficients.tsv",
        preprocessing=f"{STAGE09A_ROOT}/layer1_model_preprocessing.tsv",
        selection=f"{STAGE09A_ROOT}/layer1_model_selection.tsv",
        cv_groups=f"{STAGE09A_ROOT}/layer1_cv_by_group.tsv",
        cv_23=f"{STAGE09A_ROOT}/layer1_cv_summary_23nt.tsv",
        cv_24=f"{STAGE09A_ROOT}/layer1_cv_summary_24nt.tsv",
        representation=f"{STAGE09A_ROOT}/layer1_representation_diagnostic.tsv",
        benchmarks=f"{STAGE09A_ROOT}/layer1_architecture_benchmarks.tsv",
        candidates=f"{STAGE09A_ROOT}/candidate_layer1.tsv",
        parameters="results/09_feature_layers/stage09_parameters.tsv",
        qc="results/09_feature_layers/stage09_qc.tsv",
    params:
        legacy_core=str(LEGACY_CORE),
        candidates="results/06_targets/target_candidates.tsv",
        output_root=STAGE09A_ROOT,
    conda:
        "../envs/stage09a.yaml"
    shell:
        "python3 {input.script:q} --legacy-core {params.legacy_core:q} "
        "--candidates {params.candidates:q} --solver {input.solver:q} "
        "--model-script {input.model_script:q} --output-root {params.output_root:q}"
