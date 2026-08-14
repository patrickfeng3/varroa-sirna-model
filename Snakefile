from pathlib import Path


local_paths = Path("config/paths.local.yaml")
if local_paths.exists():
    configfile: str(local_paths)
else:
    configfile: "config/paths.example.yaml"


include: "workflow/rules/validate_legacy_core.smk"


rule all:
    input:
        "results/00_validation/legacy_core_validation.tsv",
        "results/00_validation/legacy_core_validation.md",
