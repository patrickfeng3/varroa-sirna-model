STAGE11_RESULTS = "results/11_region_explorer"
STAGE11_WEB = "web/stage11"


rule stage11_web_export:
    input:
        script="workflow/scripts/export_stage11_web_data.py",
        target_manifest="resources/targets/target_manifest.tsv",
        html=f"{STAGE11_WEB}/index.html",
        javascript=f"{STAGE11_WEB}/app.js",
        stylesheet=f"{STAGE11_WEB}/styles.css",
    output:
        web_data=f"{STAGE11_WEB}/data/Vd_CHIBIN_stage11.json",
        web_data_js=f"{STAGE11_WEB}/data/Vd_CHIBIN_stage11.js",
        qc=f"{STAGE11_RESULTS}/stage11_export_qc.tsv",
        parameters=f"{STAGE11_RESULTS}/stage11_parameters.tsv",
    params:
        stage10="results/10_candidate_integration/candidate_stage10a.tsv",
        target_id="Vd_CHIBIN",
        repo_root=".",
    shell:
        "python3 {input.script:q} --stage10 {params.stage10:q} "
        "--target-manifest {input.target_manifest:q} --target-id {params.target_id:q} "
        "--repo-root {params.repo_root:q} --web-data {output.web_data:q} "
        "--web-data-js {output.web_data_js:q} "
        "--qc {output.qc:q} --parameters {output.parameters:q}"
