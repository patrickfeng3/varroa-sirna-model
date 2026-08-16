import csv


STAGE08_ROOT = "results/08_candidate_biophysics"
STAGE08_TARGET_MANIFEST = "resources/targets/target_manifest.tsv"
STAGE08_ZUBER_RESOURCE = "resources/parameters/zuber_2022_wcf_dg37.tsv"


def stage08_transcript_fastas(_wildcards):
    with open(STAGE08_TARGET_MANIFEST, newline="", encoding="utf-8") as handle:
        return sorted({row["fasta_path"] for row in csv.DictReader(handle, delimiter="\t")})


rule stage08_candidate_biophysics:
    input:
        script="workflow/scripts/stage08.py",
        candidates="results/06_targets/target_candidates.tsv",
        target_manifest=STAGE08_TARGET_MANIFEST,
        transcript_fastas=stage08_transcript_fastas,
        zuber_resource=STAGE08_ZUBER_RESOURCE,
        analysis_config="config/analysis.yaml",
    output:
        candidates=f"{STAGE08_ROOT}/candidate_biophysics.tsv",
        parameters=f"{STAGE08_ROOT}/stage08_parameters.tsv",
        qc=f"{STAGE08_ROOT}/stage08_qc.tsv",
    params:
        output_root=STAGE08_ROOT,
    shell:
        "python3 {input.script:q} --candidates {input.candidates:q} "
        "--target-manifest {input.target_manifest:q} "
        "--zuber-resource {input.zuber_resource:q} "
        "--analysis-config {input.analysis_config:q} "
        "--output-root {params.output_root:q}"
