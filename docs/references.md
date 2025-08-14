# Overview of the reference files

References, panel of normals and design files used by the pipeline are listed below. 
They are defined in the `config/config.yaml` file and can be found in the `ref_data/` folder of the project.

| config entry / rule            | sub-entry        | file                                      | description                           |
|--------------------------------|------------------|-------------------------------------------|---------------------------------------|
| reference                      | fasta            | `GRCh38.fasta`                            | Reference genome in FASTA format      |
|                                | fai              | `GRCh38.fasta.fai`                        | Index file for the reference genome   |
|                                | trf              | `GRCh38.trf.bed`                          | Human tandem repeats                  |
|                                | design_bed       | `expected_coverage_annotated.bed`         | ???                                   |
|                                | targets_bed      | `target_bases_covered_by_probes.bed`      | ???                                   |
| cnvkit_batch                   | normal_reference | `cnvkit.PoN.cnn`                          | Panel of normals for CNVkit           |
| deep_somatic_t_only            | pon              | `noral_db.vcf.gz`                         | Panel of normals for DeepSomatic      |
| pbsv_discover                  | trf              | `GRCh38.trf.bed`                          | Human tandem repeats                  |
| hificnv                        | exclude          | `cnv.excluded_regions.common_50.hg38.bed` | regions to exclude from CNV calling   |
| bcftools_filter_include_region | panel            | `expected_coverage_annotated.bed`         | regions to include in variant calling |
| general_report                 |                  | `config/report.yaml`                      | General report configuration file     |