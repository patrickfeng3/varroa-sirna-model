STAGE03_ROOT = "results/03_steprna"


rule stage03_official_steprna:
    input:
        validator="results/00_validation/legacy_core_validation.tsv",
        script="workflow/scripts/stage03.py",
        analysis_config="config/analysis.yaml",
        environment="workflow/envs/stage03.yaml",
    output:
        accounting=f"{STAGE03_ROOT}/qc/stage03_accounting.tsv",
        software=f"{STAGE03_ROOT}/provenance/software_versions.tsv",
        run_manifest=f"{STAGE03_ROOT}/provenance/run_manifest.tsv",
        input_manifest=f"{STAGE03_ROOT}/inputs/input_manifest.tsv",
        focal_manifest=f"{STAGE03_ROOT}/inputs/focal_reference_manifest.tsv.gz",
        passenger_manifest=f"{STAGE03_ROOT}/inputs/passenger_manifest.tsv.gz",
        recovery=f"{STAGE03_ROOT}/parsed/passenger_recovery_by_pair.tsv",
        spectrum=f"{STAGE03_ROOT}/parsed/overhang_spectrum_by_pair.tsv",
        passenger_length=f"{STAGE03_ROOT}/parsed/passenger_length_by_pair.tsv",
        joint=f"{STAGE03_ROOT}/parsed/joint_geometry_by_pair.tsv",
        joint_references=f"{STAGE03_ROOT}/parsed/joint_geometry_references.tsv.gz",
    params:
        legacy_core=str(LEGACY_CORE),
        output_root=STAGE03_ROOT,
    conda:
        "../envs/stage03.yaml"
    shell:
        "python {input.script:q} --legacy-core {params.legacy_core:q} "
        "--config {input.analysis_config:q} --output-root {params.output_root:q}"


rule stage03_joint_geometry_spectrum:
    input:
        script="workflow/scripts/stage03_joint_spectrum.py",
        library="workflow/scripts/stage03.py",
        run_manifest=f"{STAGE03_ROOT}/provenance/run_manifest.tsv",
        prespecified_joint=f"{STAGE03_ROOT}/parsed/joint_geometry_by_pair.tsv",
    output:
        spectrum=f"{STAGE03_ROOT}/parsed/joint_geometry_spectrum_by_pair.tsv",
        accounting=f"{STAGE03_ROOT}/qc/stage03_joint_geometry_spectrum_accounting.tsv",
    params:
        stage03_root=STAGE03_ROOT,
    conda:
        "../envs/stage03.yaml"
    shell:
        "python {input.script:q} --stage03-root {params.stage03_root:q} "
        "--output {output.spectrum:q} --qc-output {output.accounting:q}"
