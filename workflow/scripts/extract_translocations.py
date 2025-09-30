#!/usr/bin/env python3
# this script extracts and transforms translocation data from Sniffles2 VCF files

import pandas as pd
import os
import re

# Input TSV file with sample names
input_tsv = "config/samples.tsv"
vcf_suffix = "_T.vcf"
path_prefix = "results/sniffles2_call/"
output_tsv = "results/chr4_14_breakends.tsv"
chromosome_of_interest = "chr4"

# Read sample names
samples_df = pd.read_csv(input_tsv, sep="\t")
sample_names = samples_df["sample"].tolist()

# Prepare output data
output_rows = []

for sample in sample_names:
    vcf_file = f"{path_prefix}{sample}{vcf_suffix}"
    if not os.path.exists(vcf_file):
        print(f"File not found: {vcf_file}")
        continue

    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip header lines
            parts = line.strip().split("\t")
            chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample_data = parts[:10]

            if chrom != chromosome_of_interest or "BND" not in id_:
                continue

            # --- ALT column processing ---
            # Remove square brackets
            # Remove brackets and leading/trailing base (A,C,G,T) from each allele
            alt_clean = ",".join(
                re.sub(r"^[ACGT]?([\w\d_]+:[\d]+)[ACGT]?$", r"\1", allele)
                for allele in alt.replace("[", "").replace("]", "").split(",")
            )

            # Split ALT into TRANS_CHROM and TRANS_POS
            # If ALT contains a colon, split; else leave empty
            trans_chrom, trans_pos = "", ""
            if ":" in alt_clean:
                trans_chrom, trans_pos = alt_clean.split(":", 1)

            # Parse INFO
            info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
            coverage = info_dict.get("COVERAGE", "")
            vaf = info_dict.get("VAF", "")

            # Parse FORMAT and sample data
            format_keys = format_.split(":")
            format_values = sample_data.split(":")
            format_dict = dict(zip(format_keys, format_values))

            row = {
                "SAMPLE": sample,
                "CHROM": chrom,
                "POS": pos,
                "TRANS_CHROM": trans_chrom,
                "TRANS_POS": trans_pos,
                "FILTER": filter_,
                "COVERAGE": coverage,
                "VAF": vaf,
                "GENOTYPE": format_dict.get("GT", ""),
                "GENOME QUALITY": format_dict.get("GQ", ""),
                "DEPTH REF": format_dict.get("DR", ""),
                "DEPTH TRANS": format_dict.get("DV", ""),
            }
            output_rows.append(row)

# Write to output TSV
output_df = pd.DataFrame(output_rows)
output_df.to_csv(output_tsv, sep="\t", index=False)
print(f"Output written to {output_tsv}")
