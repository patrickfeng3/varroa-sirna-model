from pathlib import Path


LEGACY_CORE = Path(config["legacy_core"]).expanduser().resolve()


rule validate_legacy_core:
    input:
        validator="workflow/scripts/validate_legacy_core.py",
    output:
        tsv="results/00_validation/legacy_core_validation.tsv",
        md="results/00_validation/legacy_core_validation.md",
    params:
        legacy_core=str(LEGACY_CORE),
    shell:
        "python3 workflow/scripts/validate_legacy_core.py "
        "--legacy-core {params.legacy_core:q} --output-tsv {output.tsv:q} --output-md {output.md:q}"
