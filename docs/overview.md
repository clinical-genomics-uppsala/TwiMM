# Overview of the pipeline
Here is a brief overview of the entire pipeline. For details see subsections and the [hydra-genetics](https://github.com/hydra-genetics/hydra-genetics) documentation.

1. **Input files**: bam, not aligned, demultiplexed 
2. **Preprocessing**: 
   - Mark duplicates with `pbmarkdup`
   - Align reads with `pbmm2`
   - Sort and index BAM files with `samtools` and `pbmm2`
3. **SNV and InDels**:
   - Call variants with `DeepSomatic`
   - Phasing and haplotagging with `whatshap`
   - Filter variants with `bcftools` and PoN
4. **CNV**:
   - Call CNVs with `cnvkit`, `pbsv`, `Sniffles2`, and `Severus`
   - Filter CNVs with `cnvkit` PoN
   - Create CNV reports with `cnvkit`