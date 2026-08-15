import csv


STAGE06_ROOT = "results/06_targets"
STAGE06_TARGET_MANIFEST = "resources/targets/target_manifest.tsv"


def stage06_registry_resources(_wildcards):
    paths = []
    with open(STAGE06_TARGET_MANIFEST, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            paths.append(row["fasta_path"])
            annotation = row["annotation_path"].strip()
            if annotation.upper() not in {"", "NA", "N/A", "NONE", "."}:
                paths.append(annotation)
    return sorted(set(paths))


rule stage06_targets:
    input:
        script="workflow/scripts/stage06.py",
        target_manifest=STAGE06_TARGET_MANIFEST,
        target_resources=stage06_registry_resources,
    output:
        reference_summary=f"{STAGE06_ROOT}/target_reference_summary.tsv",
        candidates=f"{STAGE06_ROOT}/target_candidates.tsv",
        accounting=f"{STAGE06_ROOT}/qc/stage06_accounting.tsv",
        provenance=f"{STAGE06_ROOT}/provenance/stage06_manifest.tsv",
    params:
        output_root=STAGE06_ROOT,
    shell:
        "python3 {input.script:q} --target-manifest {input.target_manifest:q} "
        "--output-root {params.output_root:q}"
