# Columns to keep and their readable names
COLUMNS_KEEP = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "Consequence", 
                       "IMPACT", "SYMBOL", "Feature_type", "BIOTYPE", "EXON", "INTRON", 
                       "Existing_variation", "VARIANT_CLASS", "AF", "GT", "GQ", "DP", "AD"]

COLUMNS_READABLE_NAMES = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "Consequence", 
                       "IMPACT", "GENE", "FEATURE", "TYPE", "EXON", "INTRON", 
                       "KNOWN_VARIATION", "VARIANT_CLASS", "AF", "GT", "GQ", "DP", "AD"]

# VEP CSQ field format
VEP_FIELDS = [
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
FORMAT_FIELDS = ["GT", "GQ", "DP", "AD", "VAF", "PL"]

# Imports
import pandas as pd
import gzip


# Functions
def parse_info(info_str):
    """
    Parse the INFO field from a VCF line
    param info_str: The INFO field from the VCF
    return: Dictionary mapping INFO keys to values
    """
    info_dict = {}
    for entry in info_str.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value
    return info_dict


def parse_format(format_str, sample_str):
    """
    Parse the FORMAT and sample fields from a VCF line
    param format_str: The FORMAT field from the VCF
    param sample_str: The sample field from the VCF
    return: Dictionary mapping FORMAT keys to sample values
    """
    keys = format_str.split(":")
    values = sample_str.split(":")
    return dict(zip(keys, values))


def parse_vcf_line(line, vep_fields, format_fields):
    """
    Parse a single VCF line into a dictionary
    param line: A line from a VCF file
    return: Dictionary with parsed fields
    """
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
    """
    Open a VCF file, handling gzip if necessary
    param vcf_path: Path to the VCF file
    return: File object
    """
    if vcf_path.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    else:
        return open(vcf_path, "r")


def vcf_to_df(vcf_path, vep, fields):
    """
    Convert a VEP annotated VCF file to a DataFrame
    param vcf_path: Path to the VCF file
    return: DataFrame with parsed VCF data
    """
    rows = []
    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            row = parse_vcf_line(line, vep, fields)
            rows.append(row)
    df = pd.DataFrame(rows)
    return df
    

def pick_vcf_columns(vcf_df: pd.DataFrame, columns_to_keep: list = None) -> pd.DataFrame:
    """
    Pick relevant columns from the VCF DataFrame
    param vcf_df: DataFrame with VCF data
    return: DataFrame with selected columns
    """   
    return vcf_df[columns_to_keep]


def filter_vcf(vcf_df, column, value):
    """
    Filter VCF DataFrame based on value in column.
    param vcf_df: DataFrame with VCF data
    param column: Column name to filter on
    param value: Value to filter by
    return: Filtered DataFrame
    """
    return vcf_df[vcf_df[column] == value]


if __name__ == "__main__":
    # Get input and output paths from snakemake
    vcf_file=snakemake.input.vcf
    output_xlsx=snakemake.output.xlsx

    vcf_df=pick_vcf_columns(vcf_to_df(vcf_file, VEP_FIELDS, FORMAT_FIELDS), COLUMNS_KEEP)
    vcf_df.columns=COLUMNS_READABLE_NAMES

    snv_tp53=filter_vcf(vcf_df, "GENE", "TP53")

    snv_df=filter_vcf(vcf_df, "FILTER", "PASS")

    with pd.ExcelWriter(output_xlsx, engine='xlsxwriter') as writer:
        # All the SNVs in one sheet
        snv_df.to_excel(writer, sheet_name='SNV', index=False)
        workbook  = writer.book
        worksheet_snv = writer.sheets['SNV']
        # Get the dimensions of the dataframe
        max_row, max_col = snv_df.shape
        # Set autofilter on the header row
        worksheet_snv.autofilter(0, 0, max_row, max_col - 1)
        
        # TP53 SNVs in a separate sheet
        snv_tp53.to_excel(writer, sheet_name='TP53', index=False)
        # separate sheet for TP53 variantsß
        worksheet_tp53 = writer.sheets['TP53']
        # its dimensions
        max_row, max_col = snv_tp53.shape
        # Set autofilter on the header row 
        worksheet_tp53.autofilter(0, 0, max_row, max_col - 1)
