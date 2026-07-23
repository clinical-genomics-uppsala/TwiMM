# Overview of the reference files

References, panel of normals and design files used by the pipeline are listed below. 
They are defined in the `config/config.yaml` file and can be found in the `ref_data/` folder of the project.

| config entry / rule            | sub-entry        | file                                      | description                                                                       |
|--------------------------------|------------------|-------------------------------------------|-----------------------------------------------------------------------------------|
| reference                      | fasta            | `GRCh38.fasta`                            | Reference genome in FASTA format                                                  |
|                                | fai              | `GRCh38.fasta.fai`                        | Index file for the reference genome                                               |
|                                | trf              | `GRCh38.trf.bed`                          | Human tandem repeats (used by Severus and PBSV)                                   |
|                                | design_bed       | `expected_coverage_annotated.bed`         | BED file with panel target regions                                                |
|                                | sv_databases     | `gnomad.v4.1.sv.sites.no_cnv.PASS.vcf.gz`, `sv_normal.vcf.gz` | Population SV VCF files queried directly by SVDB. List — add new population databases here. |
| cnvkit_batch                   | normal_reference | `cnvkit.PoN.cnn`                          | Panel of normals for CNVkit                                                       |
| deepsomatic_t_only             | pon              | `snv_normal.vcf.gz`                       | Panel of normals for DeepSomatic (used when `use_deepsomatic: true`)              |
| severus_t_only                 | pon              | `PoN_1000G_hg38.tsv.gz`                   | Panel of normals for Severus (1000 Genomes project)                               |
|                                | vntr             | `GRCh38.trf.bed`                          | Tandem repeat BED file passed to Severus via `--vntr`                             |
| bcftools_filter_include_region | panel            | `expected_coverage_annotated.bed`         | Regions to include in variant calling                                             |
| general_report                 |                  | `config/report.yaml`                      | General report configuration file                                                 |
| merge_cnv_json                 | cytobands        | `config/cytoBand.hg38.txt`                | UCSC cytoband definitions (hg38) for the CNV HTML report's chromosome plot        |
| merge_cnv_json                 | table_filter_config | `config/filters/table_filter.yaml`     | CNV classification/filter rules for the CNV HTML report's results table           |