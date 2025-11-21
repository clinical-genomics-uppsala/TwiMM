# Imports
import pandas as pd
import gzip
import re
import logging
import yaml


# Functions
def parse_info(info_str: str) -> dict:
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


def parse_format(format_str: str, sample_str: str) -> dict:
    """
    Parse the FORMAT and sample fields from a VCF line
    param format_str: The FORMAT field from the VCF
    param sample_str: The sample field from the VCF
    return: Dictionary mapping FORMAT keys to sample values
    """
    keys = format_str.split(":")
    values = sample_str.split(":")
    return dict(zip(keys, values))


def parse_vcf_line(line, vep_fields: list, format_fields: list) -> dict:
    """
    Parse a single VCF line into a dictionary
    param line: A line from a VCF file
    param vep_fields: List of VEP annotation fields
    param format_fields: List of FORMAT fields to extract
    return: Dictionary with parsed fields
    """
    fields = line.strip().split("\t")
    chrom, pos, _, ref, alt, qual, fltr, info, fmt, sample = fields[:10]

    info_dict = parse_info(info)
    format_dict = parse_format(fmt, sample)

    csq_data = info_dict.get("CSQ", "").split(",")[0]  # take first annotation
    csq_values = csq_data.split("|")
    csq_dict = dict(zip(vep_fields, csq_values))

    row = {
        "CHROM": chrom,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "QUAL": qual,
        "FILTER": fltr,
    }

    row.update(csq_dict)
    for key in format_fields:
        row[key] = format_dict.get(key, "")

    return row


def parse_sv_vcf_line(line: str) -> dict:
    """
    Parse a single structural variant VCF line into a dictionary
    param line: A line from a VCF file
    return: Dictionary with parsed fields
    """
    parts = line.strip().split("\t")
    chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample_data = parts[:10]

    # Clean ID field
    match = re.search(r"\.(.*?)\.", id_)
    id_clean = match.group(1)

    # Parse INFO
    info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
    coverage = info_dict.get("COVERAGE", "")
    vaf = info_dict.get("VAF", "")
    support = info_dict.get("SUPPORT", "")
    strand = info_dict.get("STRAND", "")
    # available for INS & DEL only
    end = info_dict.get("END", "")
    svlen = info_dict.get("SVLEN", "")

    # Parse FORMAT and sample data
    format_keys = format_.split(":")
    format_values = sample_data.split(":")
    format_dict = dict(zip(format_keys, format_values))

    # this may require subsetting depending on your needs
    row = {
        "CHROM": chrom,
        "POS": pos,
        "END": end,
        "TYPE": id_clean,
        "SVLEN": svlen,
        "ALT": alt,
        "FILTER": filter_,
        "COVERAGE": coverage,
        "SUPPORT": support,
        "STRAND": strand,
        "VAF": vaf,
        "GENOTYPE": format_dict.get("GT", ""),
        "GENOME QUALITY": format_dict.get("GQ", ""),
        "DEPTH REF": format_dict.get("DR", ""),
        "DEPTH TRANS": format_dict.get("DV", ""),
    }
    return row


def parse_cnvkit_vcf_line(vcf_line: str) -> dict:
    """
    Parse a single CNVkit VCF line into a dictionary
    param vcf_line: A line from a CNVkit VCF file
    return: Dictionary with parsed fields
    """
    # Split the line into its tab-separated fields
    fields = vcf_line.strip().split("\t")

    # Basic columns
    chrom = fields[0]
    pos = int(fields[1])
    variant_type = fields[4].strip("<>")  # Remove angle brackets from ALT field

    # Parse INFO field
    info_dict = parse_info(fields[7])

    # Extract desired INFO fields
    genes = info_dict.get("Genes", "")
    end = int(info_dict.get("END", 0))
    svlen = int(info_dict.get("SVLEN", 0))
    log_odds_ratio = float(info_dict.get("LOG_ODDS_RATIO", "nan"))
    corr_cn = float(info_dict.get("CORR_CN", "nan"))
    probes = int(info_dict.get("PROBES", 0))
    try:
        baf = float(info_dict.get("BAF", "nan"))
    except ValueError:
        logging.info(
            f"Non-numeric BAF value found: {info_dict.get('BAF')}, will use it as-is"
        )
        baf = info_dict.get("BAF", "nan")

    # Parse FORMAT and sample fields
    format_dict = parse_format(format_str=fields[8], sample_str=fields[9])
    # extract certain fields
    gt = format_dict.get("GT", "")
    cnq = float(format_dict.get("CNQ", ""))
    dp = float(format_dict.get("DP", ""))

    # Build the row dictionary
    row = {
        "CHROM": chrom,
        "POS": pos,
        "VARIANT_TYPE": variant_type,
        "GENE": genes,
        "END": end,
        "SVLEN": svlen,
        "LOG_ODDS_RATIO": log_odds_ratio,
        "CORR_CN": corr_cn,
        "PROBES": probes,
        "BAF": baf,
        "GT": gt,
        "CNQ": cnq,
        "DP": dp,
    }

    return row


def open_vcf(vcf_path: str):
    """
    Open a VCF file, handling gzip if necessary
    param vcf_path: Path to the VCF file
    return: File object
    """
    if vcf_path.endswith(".gz"):
        return gzip.open(vcf_path, "rt")
    else:
        return open(vcf_path, "r")


def vcf_to_df(vcf_path: str, vep: list, fields: list) -> pd.DataFrame:
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


def sv_vcf_to_df(vcf_path: str, cnvkit: bool) -> pd.DataFrame:
    """
    Convert a structural variant VCF file to a DataFrame
    param vcf_path: Path to the VCF file
    param cnvkit: whether to parse CNVkit-specific fields
    return: DataFrame with parsed VCF data
    """
    parse_line = parse_cnvkit_vcf_line if cnvkit else parse_sv_vcf_line

    with open_vcf(vcf_path) as vcf:
        rows = [parse_line(line) for line in vcf if not line.startswith("#")]

    return pd.DataFrame(rows)


def pick_vcf_columns(
    vcf_df: pd.DataFrame, columns_to_keep: list = None
) -> pd.DataFrame:
    """
    Pick relevant columns from the VCF DataFrame
    param vcf_df: DataFrame with VCF data
    return: DataFrame with selected columns
    """
    return vcf_df[columns_to_keep]


if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(
        filename=snakemake.log[0],
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
        level=logging.INFO,
    )

    logging.info("Script started")
    logging.info(f"Sample name: {snakemake.wildcards.sample}")

    # Get input and output paths from snakemake
    vcf_snv = snakemake.input.vcf_snv
    vcf_sv = snakemake.input.vcf_sv
    vcf_cnv = snakemake.input.vcf_cnv
    output_xlsx = snakemake.output.xlsx

    logging.info(
        f"Input files: SNV VCF: {vcf_snv}, SV VCF: {vcf_sv}, CNV VCF: {vcf_cnv}\nOutput file: {output_xlsx}"
    )

    # get params as lists
    filter_yaml_file = snakemake.params.filter_config
    with open(filter_yaml_file) as file:
        filters = yaml.load(file, Loader=yaml.FullLoader)
    format_fields = filters.get("format_fields", [])
    vep_fields = filters.get("vep_info_fields", [])
    columns_keep = filters.get("columns_keep", [])
    snvs_remove = filters.get("snvs_remove", [])
    idid_min_len = filters.get("idid_min_len", 1000)

    if any(
        x is None
        for x in [
            snvs_remove,
            format_fields,
            vep_fields,
            columns_keep,
            idid_min_len,
        ]
    ):
        logging.error("Missing parameters")
        raise ValueError(
            "Some required parameters are missing. Check your config file."
        )

    # read SNV vcf file
    logging.info("Reading provided VCF files")
    snv_all_df = vcf_to_df(vcf_snv, vep_fields, format_fields)
    
    # remove not important SNV categories and those not passing default filter
    snv_all_df = snv_all_df[
        (~snv_all_df["Consequence"].isin(snvs_remove)) & (snv_all_df["FILTER"] == "PASS")
    ]

    # keep only chosen columns
    snv_picked_columns = pick_vcf_columns(snv_all_df, columns_keep)

    # rename SYMBOL to GENE for clarity
    snv_picked_columns = snv_picked_columns.rename(columns={"SYMBOL": "GENE"})
    
    # Collect TP53 SNV to a separate dataframe
    snv_tp53 = snv_picked_columns[snv_picked_columns["GENE"] == "TP53"]
    logging.info(f"TP53 SNVs after filtering: {len(snv_tp53)}")

    # Collect the rest of SNVs to a separate dataframe
    snv_rest = snv_picked_columns[snv_picked_columns["GENE"] != "TP53"]
    logging.info(f"Not TP53 SNVs after filtering: {len(snv_rest)}")

    # read SV vcf file
    sv_df = sv_vcf_to_df(vcf_sv, cnvkit=False)
    logging.info(f"Total SVs read: {len(sv_df)}")
    
    # filter both chr4 and BND
    tn_chr4 = sv_df[(sv_df["CHROM"] == "chr4") & (sv_df["TYPE"] == "BND")]
    logging.info(f"Translocations from chr4: {len(tn_chr4)}")
    
    # filter both chr14 and BND
    tn_chr14 = sv_df[(sv_df["CHROM"] == "chr14") & (sv_df["TYPE"] == "BND")]
    logging.info(f"Translocations from chr14: {len(tn_chr14)}")

    # read SV vcf file and extract IDID variants on chr14
    sv_chr14_pass = sv_df[(sv_df["CHROM"] == "chr14") & (sv_df["FILTER"] == "PASS")]
    # convert SVLEN to numeric and turn empty strings to NaN
    sv_chr14_pass["SVLEN"] = pd.to_numeric(sv_chr14_pass["SVLEN"], errors="coerce")
    # keep TYPE!=BND
    sv_chr14_idid = sv_chr14_pass[
        (~sv_chr14_pass["TYPE"].isin(["BND"]))
        & (sv_chr14_pass["SVLEN"].abs() >= idid_min_len)
    ]
    logging.info(f"Total IDID variants on chr14: {len(sv_chr14_idid)}")

    # read CNVkit VCF file
    cnv_df = sv_vcf_to_df(vcf_cnv, cnvkit=True)
    logging.info(f"Total CNVs read: {len(cnv_df)}")

    with pd.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        # All the SNVs in one sheet
        snv_rest.to_excel(writer, sheet_name="SNVs", index=False)
        workbook = writer.book
        worksheet_snv = writer.sheets["SNVs"]
        # Get the dimensions of the dataframe
        max_row, max_col = snv_rest.shape
        # Set autofilter on the header row
        worksheet_snv.autofilter(0, 0, max_row, max_col - 1)

        # TP53 SNVs in a separate sheet
        snv_tp53.to_excel(writer, sheet_name="TP53", index=False)
        # separate sheet for TP53 variantsß
        worksheet_tp53 = writer.sheets["TP53"]
        # its dimensions
        max_row, max_col = snv_tp53.shape
        # Set autofilter on the header row
        worksheet_tp53.autofilter(0, 0, max_row, max_col - 1)

        # Translocations (BND) from chr4
        tn_chr4.to_excel(writer, sheet_name="Tn_chr4", index=False)
        worksheet_sv = writer.sheets["Tn_chr4"]
        max_row, max_col = tn_chr4.shape
        worksheet_sv.autofilter(0, 0, max_row, max_col - 1)

        # Translocations (BND) from chr14
        tn_chr14.to_excel(writer, sheet_name="Tn_chr14", index=False)
        worksheet_sv = writer.sheets["Tn_chr14"]
        max_row, max_col = tn_chr14.shape
        worksheet_sv.autofilter(0, 0, max_row, max_col - 1)

        # IDID variants from chr14
        sv_chr14_idid.to_excel(writer, sheet_name="IDID_chr14", index=False)
        worksheet_idid = writer.sheets["IDID_chr14"]
        max_row, max_col = sv_chr14_idid.shape
        worksheet_idid.autofilter(0, 0, max_row, max_col - 1)

        # CNVs in a separate sheet
        cnv_df.to_excel(writer, sheet_name="CNV", index=False)
        worksheet_cnv = writer.sheets["CNV"]
        max_row, max_col = cnv_df.shape
        worksheet_cnv.autofilter(0, 0, max_row, max_col - 1)

    logging.info("Script finished successfully")
