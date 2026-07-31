# Software used in the TwiMM pipeline

The TwiMM pipeline utilizes various software tools to process and analyze PacBio hybrid capture sequencing data. 
Below is a list of the key software components used in the pipeline:

| Software       | Version    | Description                                                                                                                  |
|----------------|------------|------------------------------------------------------------------------------------------------------------------------------|
| snakemake      | >=7.32.4,<8 | A workflow management system that allows for the creation of reproducible and scalable data analyses.                        |
| hydra-genetics | 4.0.0      | A framework for building and managing Snakemake pipelines, providing a structured approach to bioinformatics workflows.      |
| DeepSomatic    | 1.9.0      | A tool for somatic variant calling in long-read sequencing data, specifically designed for PacBio data. Used when `use_deepsomatic: true`. |
| ClairS-TO      | 0.4.4      | Tumor-only somatic small variant caller for long-read data. Default SNV caller (`use_deepsomatic: false`).                   |
| bcftools       | 1.16       | A set of utilities for manipulating variant call format (VCF) and binary call format (BCF) files.                            |
| samtools       | 1.16       | A suite of programs for interacting with high-throughput sequencing data, including BAM file manipulation.                   |
| bedtools       | 2.30.0     | A powerful toolset for genome arithmetic, allowing for the manipulation of genomic intervals and annotations.                |
| CNVkit         | 0.9.11     | A toolkit for detecting copy number variations in genomic data, particularly useful for analyzing long-read sequencing data. |
| Severus        | 1.5        | Tumor-only SV caller optimized for long reads, supporting haplotagged BAMs and panel of normals.                             |
| PBSV           | 2.9.0      | PacBio SV caller using split reads and soft clips from pbmm2-aligned BAMs.                                                  |
| Sniffles2      | 2.6.1      | A structural variant caller for long-read sequencing data, providing accurate detection of large genomic rearrangements.     |
| SVDB           | 2.8.4      | SV database tool used to build a population frequency database (`--build`) and annotate sample SVs (`--query`).             |
| whatshap       | 2.2        | A tool for phasing and haplotagging variants in sequencing data.                                                             |
| pbmarkdup      | 1.1.0      | A tool for marking duplicates in PacBio sequencing data, helping to improve the accuracy of variant calling.                 |
| pbmm2          | 1.16       | A fast and efficient aligner for PacBio sequencing data, providing high-quality alignments for downstream analysis.          |
| VACmap         | 1.0.0      | An alternative long-read aligner; selected when `aligner: "vacmap"` is set in config.                                       |
