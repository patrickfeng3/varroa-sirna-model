STAGE04_ROOT = "results/04_duplex_geometry"


rule stage04_duplex_geometry:
    input:
        script="workflow/scripts/stage04.py",
        analysis_config="config/stage04.yaml",
        stage02_expected="results/02_terminal_enrichment/background/terminal_expected_by_pair.tsv",
        stage02_pair="results/02_terminal_enrichment/enrichment/terminal_enrichment_by_pair.tsv",
        stage02_across="results/02_terminal_enrichment/enrichment/terminal_enrichment_across_dataset.tsv",
        stage03_runs="results/03_steprna/provenance/run_manifest.tsv",
        stage03_focals="results/03_steprna/inputs/focal_reference_manifest.tsv.gz",
        stage03_recovery="results/03_steprna/parsed/passenger_recovery_by_pair.tsv",
        stage03_spectrum="results/03_steprna/parsed/overhang_spectrum_by_pair.tsv",
        stage03_joint="results/03_steprna/parsed/joint_geometry_by_pair.tsv",
        stage03_joint_refs="results/03_steprna/parsed/joint_geometry_references.tsv.gz",
    output:
        accounting=f"{STAGE04_ROOT}/qc/stage04_accounting.tsv",
        full_sample=f"{STAGE04_ROOT}/population/full_spectrum_by_sample.tsv",
        full_across=f"{STAGE04_ROOT}/population/full_spectrum_across_dataset.tsv",
        recovery=f"{STAGE04_ROOT}/population/passenger_recovery_across_dataset.tsv",
        joint_sample=f"{STAGE04_ROOT}/population/joint_geometry_by_sample.tsv",
        joint_across=f"{STAGE04_ROOT}/population/joint_geometry_across_dataset.tsv",
        paired=f"{STAGE04_ROOT}/comparisons/paired_23_vs_24.tsv",
        sequence_pair=f"{STAGE04_ROOT}/sequence_features/geometry_terminal_by_pair.tsv",
        sequence_sample=f"{STAGE04_ROOT}/sequence_features/geometry_terminal_by_sample.tsv",
        sequence_across=f"{STAGE04_ROOT}/sequence_features/geometry_terminal_across_dataset.tsv",
        contrasts=f"{STAGE04_ROOT}/sequence_features/geometry_specific_contrasts.tsv",
        redundancy=f"{STAGE04_ROOT}/sequence_features/redundancy.tsv",
    params:
        stage02_root="results/02_terminal_enrichment",
        stage03_root="results/03_steprna",
        output_root=STAGE04_ROOT,
    shell:
        "python3 {input.script:q} --stage02-root {params.stage02_root:q} "
        "--stage03-root {params.stage03_root:q} --config {input.analysis_config:q} "
        "--output-root {params.output_root:q}"


rule stage04_joint_geometry_spectrum:
    input:
        script="workflow/scripts/stage04_joint_spectrum.py",
        library="workflow/scripts/stage04.py",
        analysis_config="config/stage04.yaml",
        pair_spectrum="results/03_steprna/parsed/joint_geometry_spectrum_by_pair.tsv",
    output:
        sample=f"{STAGE04_ROOT}/population/joint_geometry_spectrum_by_sample.tsv",
        across=f"{STAGE04_ROOT}/population/joint_geometry_spectrum_across_dataset.tsv",
        modes=f"{STAGE04_ROOT}/population/joint_geometry_mode_by_pair.tsv",
        summary=f"{STAGE04_ROOT}/population/joint_geometry_spectrum_summary.tsv",
        accounting=f"{STAGE04_ROOT}/qc/stage04_joint_geometry_spectrum_accounting.tsv",
    params:
        output_root=STAGE04_ROOT,
    shell:
        "python3 {input.script:q} --stage03-spectrum {input.pair_spectrum:q} "
        "--config {input.analysis_config:q} --output-root {params.output_root:q}"
