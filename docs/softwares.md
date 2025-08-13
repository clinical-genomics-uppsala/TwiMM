# Software used in the TwiMy pipeline

The TwiMy pipeline utilizes various software tools to process and analyze PacBio hybrid capture sequencing data. 
Below is a list of the key software components used in the pipeline:

| Software       | Version    | Description                                                                                                                  |
|----------------|------------|------------------------------------------------------------------------------------------------------------------------------|
| snakemake      | >=7.8.0,<8 | A workflow management system that allows for the creation of reproducible and scalable data analyses.                        |
| hydra-genetics | 3.0.0      | A framework for building and managing Snakemake pipelines, providing a structured approach to bioinformatics workflows.      |
| deepsomatic    | 1.8.0      | A tool for somatic variant calling in long-read sequencing data, specifically designed for PacBio data.                      |
| bcftools       | 1.16       | A set of utilities for manipulating variant call format (VCF) and binary call format (BCF) files.                            |
| samtools       | 1.16       | A suite of programs for interacting with high-throughput sequencing data, including BAM file manipulation.                   |
| bedtools       | 2.30.0     | A powerful toolset for genome arithmetic, allowing for the manipulation of genomic intervals and annotations.                |
| CNVkit         | 0.9.9      | A toolkit for detecting copy number variations in genomic data, particularly useful for analyzing long-read sequencing data. |
| HiFiCNV        | 1.0.1      | A tool specifically designed for detecting copy number variations in HiFi (High Fidelity) long-read sequencing data.         |
| Sniffles2      | 2.6.1      | A structural variant caller for long-read sequencing data, providing accurate detection of large genomic rearrangements.     |
| pbsv           | 2.9.0      | A tool for processing and analyzing PacBio sequencing data, including variant calling and structural variant detection.      |
| whatshap       | 2.2        | A tool for inferring haplotypes from sequencing data, useful for understanding genetic variation in populations.             |
| pbmarkdup      | 1.1.0      | A tool for marking duplicates in PacBio sequencing data, helping to improve the accuracy of variant calling.                 |
| pbmm2          | 1.16       | A fast and efficient aligner for PacBio sequencing data, providing high-quality alignments for downstream analysis.          |
