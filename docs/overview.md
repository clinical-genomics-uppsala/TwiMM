# Overview of the pipeline

Here is a brief overview of the entire pipeline. For details see subsections and the [hydra-genetics](https://github.com/hydra-genetics/hydra-genetics) documentation.

1. **Input files**: bam, not aligned, demultiplexed  


2. **Preprocessing**:   
  
   2.1 Mark duplicates with `pbmarkdup`.
  
   2.2 Align reads with `pbmm2` or `VACmap` (configured via `aligner` key).

3. **SNV and InDels**:  
  
   3.1 Call variants with `DeepSomatic` (tumor-only mode) or `ClairS-TO`. 
  
   3.2 Phasing with `whatshap` (uses `whatshap phase`).
  
   3.3 Haplotagging with `whatshap` (uses `whatshap haplotag`) - produces haplotagged BAMs.
  
   3.4 Annotation of variants using `VEP`.
  
   3.5 Filter variants with `bcftools` (based on "germline" checks from panel of normals).

4. **Structural Variants (SV)**:

   4.1 Call SVs in parallel with three callers: `Severus` (tumor-only, uses haplotagged BAM and panel of normals), `PBSV` (uses pbmm2-aligned BAM), and `Sniffles2` (uses haplotagged BAM).

   4.2 Filter SV calls per caller to panel regions with `bcftools`.

   4.3 Build a population SV database (once globally) from population VCF files (gnomAD SV, custom PoN) using `svdb --build`.

   4.4 Merge per-sample SV calls from all three callers into a single VCF using `svdb --merge` (priority: Severus > Sniffles2 > PBSV).

   4.5 Annotate the merged VCF with population SV frequencies using `svdb --query` against the pre-built population database.

5. **Copy Number Variants (CNV)**:
   
   5.1 Call CNVs with `cnvkit` (using haplotagged BAMs).
   
   5.2 Annotate CNVs with `annotate_cnv` (identifying specific gene overlaps).

6. **QC and Depth of Coverage**:
   
   6.1 Calculate depth of coverage with `mosdepth`.

7. **Reporting**:
   
   7.1 Create HTML reports with `cnvkit`.
   
   7.2 Create Excel reports with combined data on SNV, CNV, and SV, including a **Software Versions** tab listing the tool versions used.
