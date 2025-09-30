import csv
import os
import pandas as pd
from tqdm import tqdm
import gzip


# VEP CSQ field format
vep_fields = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "VARIANT_CLASS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "AF",
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_MID_AF",
    "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
]

# FORMAT fields to extract
format_fields = ["GT", "GQ", "DP", "AD", "VAF", "PL"]


def parse_info(info_str):
    info_dict = {}
    for entry in info_str.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value
    return info_dict


def parse_format(format_str, sample_str):
    keys = format_str.split(":")
    values = sample_str.split(":")
    return dict(zip(keys, values))


def parse_vcf_line(line):
    fields = line.strip().split("\t")
    chrom, pos, _, ref, alt, qual, fltr, info, fmt, sample = fields[:10]

    info_dict = parse_info(info)
    format_dict = parse_format(fmt, sample)

    csq_data = info_dict.get("CSQ", "").split(",")[0]  # take first annotation
    csq_values = csq_data.split("|")
    csq_dict = dict(zip(vep_fields, csq_values))

    row = {"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, "QUAL": qual, "FILTER": fltr}

    row.update(csq_dict)
    for key in format_fields:
        row[key] = format_dict.get(key, "")

    return row


def open_vcf(vcf_path):
    if vcf_path.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    else:
        return open(vcf_path, "r")


def vcf_to_excel(vcf_path, excel_path):
    rows = []
    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            row = parse_vcf_line(line)
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_excel(excel_path, index=False)


# Input TSV file with sample names
input_samples = "config/samples.tsv"
# vcf_suffix = "_T.deepsomatic.phased.vep_annotated.vcf"
vcf_suffix = "_T.clairs.phased.include.panel.vep_annotated.vcf.gz"
path_prefix = "results/clairs_to/"
output_dir = "results/clairs_to/"

# Read sample names
samples_df = pd.read_csv(input_samples, sep="\t")
sample_names = samples_df["sample"].tolist()

for sample in tqdm(sample_names):
    vcf_file = f"{path_prefix}{sample}{vcf_suffix}"
    if not os.path.exists(vcf_file):
        print(f"File not found: {vcf_file}")
        continue

    vcf_to_excel(vcf_file, f"{output_dir}{sample}_T.clairs.phased.include.panel.vep_annotated.xlsx")
