---
name: TwiMM
description: "TwiMM is a Snakemake-based bioinformatics workflow. Its primary purpose is to analyse hybrid capture long-read (PacBio HiFi) sequencing data from the multiple myeloma gene panel. It detects SNVs/InDels, structural variants (SV), and copy number variants (CNV). SV calling uses three callers in parallel (Severus, PBSV, Sniffles2) whose outputs are merged and annotated with population frequencies via SVDB."
---

## Overview

**TwiMM** is a Snakemake-based bioinformatics pipeline for analysing hybrid-capture long-read (PacBio HiFi) sequencing data from a multiple myeloma gene panel. It is built on the [hydra-genetics](https://github.com/hydra-genetics/) framework, which provides modular, reusable rule collections loaded dynamically from GitHub at runtime.

## Environment & Package Management

All commands are run via **pixi** (the conda-based environment manager defined in `pixi.toml`). Do not use pip or conda directly.

Key tasks:

| Task | Command |
|---|---|
| Dry run | `pixi run all-dry` |
| Full run (SLURM) | `pixi run all-full` |
| Test dry run | `pixi run test-dry` |
| Test full run (SLURM) | `pixi run test-full` |
| Unit tests | `pixi run tests` |
| Pre-commit checks | `pixi run precom` |
| Format Snakefiles | `pixi run fmt-wf` |

## Architecture

### Entry points

- `workflow/Snakefile` – imports `rules/common.smk`, declares the `all` rule, then imports and wires up hydra-genetics modules using Snakemake's `module` + `use rule … as …` syntax.
- `workflow/rules/common.smk` – loads and validates all config/samples/units files, defines wildcard constraints, and provides helper functions (`get_tc_file`, `get_snv_caller_output`, `get_local_vcfs_for_svdb_merge`, etc.).

### Hydra-genetics modules

All Hydra-genetics repositories can be found here: https://github.com/orgs/hydra-genetics/repositories

The pipeline pulls modular rule sets from GitHub at runtime (versions pinned in `config/config.yaml` under `module_versions`):

| Module | Purpose |
|---|---|
| `prealignment` | Mark duplicates in CCS reads (pbmarkdup) |
| `alignment` | Map reads: pbmm2 or VACmap; merge & index BAMs |
| `snv_indels` | SNV/indel calling (ClairS-TO **or** DeepSomatic-TO), WhatsHap phase+haplotag |
| `cnv_sv` | CNVkit (CNV), PBSV + Sniffles2 + Severus (SV), SVDB merge/query |
| `annotation` | VEP annotation, `annotate_cnv` |
| `filtering` | bcftools region filter, `filter_vcf` |
| `qc` | mosdepth, sequali, Picard metrics, MultiQC |
| `reports` | CNV HTML report (JSON → HTML via D3 templates in `workflow/templates/`) |

Modules are stored in snakemake cache here: /Users/andgu885/Library/Caches/snakemake/snakemake/source-cache/
The cache has to be cleaned if the modules have been updated on the remote and this causes issues during the pipeline run.

### Configuration system

The pipeline requires at least one `--configfile`. Standard usage:

- **Full run**: `config/config.yaml` only
- **Test run**: `config/config.yaml` + `config/config_test.yaml`

Key top-level config options:

- `use_deepsomatic: true/false` – switches the SNV caller between ClairS-TO (default) and DeepSomatic-TO
- `aligner: "vacmap"` – selects the aligner (vacmap or pbmm2)
- `tc_method` (within `svdb_merge`) – defines tumor-content method; currently `pathology` (purity from `samples.tsv`)
- `reference.sv_databases` – list of population SV VCF files (e.g. gnomAD SV, custom PoN) queried directly by SVDB

### SV calling and SVDB

SV calling uses three callers in parallel: **Severus**, **PBSV**, and **Sniffles2** (all run on haplotagged BAMs). Their outputs are combined and annotated through two SVDB steps:

1. **`cnv_sv_svdb_merge`** (from `cnv_sv` module) — per-sample. Merges VCFs from all three SV callers into a single VCF using caller priority (`severus > sniffles2 > pbsv`).
2. **`cnv_sv_svdb_query`** (from `cnv_sv` module) — per-sample. Annotates the merged caller VCF with SV frequencies queried directly from the VCF files listed in `reference.sv_databases` (passed as `--db`). Querying VCF files directly (rather than a pre-built SQLite database) preserves per-source allele frequency fields, allowing gnomAD and custom PoN frequencies to be distinguished.

To add a new population SV database, append its path to `reference.sv_databases` in `config/config.yaml`.

### Input files

- `config/samples.tsv` – one row per sample; columns: `sample` (+ optional `TC` for tumor content)
- `config/units.tsv` – one row per BAM file; columns: `sample`, `type` (T/N), `platform`, `machine`, `processing_unit`, `barcode`, `methylation`, `bam`

This pipeline is designed for PacBio HiFi reads only.
Units for PacBio are indexed by `(sample, type, processing_unit, barcode)`.

Do not run the pipeline on real data on this computer. Only test dry-run is allowed.

### Output specification

`config/output_files.yaml` declares final outputs copied to `results/`. The schema is validated against `workflow/schemas/output_files.schema.yaml`.

### Local scripts

- `workflow/scripts/compile_xlsx_report.py` – combines SNV, SV, and CNV VCFs into a per-sample Excel report; invoked as a Snakemake `script:` rule. Includes a **Software Versions** tab populated from container image tags defined in `config.yaml` (aligners, SV/CNV callers, and the active SNV caller).

### Containers
Container images for this pipeline can be found at https://hub.docker.com/u/hydragenetics

## Code Style

- Snakemake files: formatted with **snakefmt** (line length 130), config in `.github/linters/.snakefmt.toml`
- Python scripts: **pycodestyle** with `--max-line-length=130`
- Commits must follow the **Conventional Commits** specification (enforced by `conventional-prs.yml`)
- Main development branch is `develop`

## Development
Add new rules from `hydra-genetics` repos. If it's not possible, write local rules.
After adding a new rule, update schemas and notify the user that rulegraphs can only be updated on the server.
