# Overview of the pipeline

Here is a brief overview of the entire pipeline. For details see subsections and the [hydra-genetics](https://github.com/hydra-genetics/hydra-genetics) documentation.

1. **Input files**: bam, not aligned, demultiplexed  
2. **Preprocessing**:   
  2.1 Mark duplicates with `pbmarkdup`.
  2.2 Align reads with `pbmm2`.
  2.3 Sort and index BAM files with `samtools` and `pbmm2`.
3. **SNV and InDels**:  
  3.1 Call variants with `DeepSomatic` (tumor-only mode) or `ClairS-TO`.
   3.2 Phasing and haplotagging with `whatshap`.
   3.3 Filter variants with `bcftools` and custom panel of normals.
   3.4 Annotation of variants using `VEP`.
4. **CNV**:
   4.1 Call CNVs with `cnvkit`, and `Sniffles2`.
   4.2 Filter CNVs with `cnvkit` PoN.
5. **Reporting**:
   5.1 Create HTML reports with `cnvkit`.
   5.2 Create Excel reports with combined data on SNV, CNV, and SV.
