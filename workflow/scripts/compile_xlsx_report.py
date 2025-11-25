# Imports
import pandas as pd
import gzip
import re
import logging
import yaml
from typing import Callable, TextIO
from pathlib import Path
from contextlib import contextmanager


# Functions
def parse_info(info_str: str) -> dict[str, str]:
    """
    Parse the INFO field from a VCF line

    Args:
        info_str: The INFO field from the VCF, e.g. "FAU=46;FCU=28"

    Returns:
        Dictionary mapping INFO keys to values, e.g. {"FAU": "46", "FCU": "28"}
    """
    return dict(entry.split("=", 1) for entry in info_str.split(";") if "=" in entry)


def validate_info_values(info_dict: dict[str, str]) -> dict[str, int]:
    """
    Validate that all values in the INFO dictionary are integers.

    Args:
        info_dict: Dictionary from parse_info.

    Returns:
        Dictionary with integer values.

    Raises:
        ValueError: If any value is not an integer.
    """
    validated = {}
    for key, value in info_dict.items():
        if not value.isdigit():
            raise ValueError(f"Invalid value for {key}: '{value}' (must be integer)")
        validated[key] = int(value)
    return validated


def validate_info_structure(info_str: str) -> None:
    """
    Validate that all entries in the INFO string are proper key=value pairs.

    Args:
        info_str: Raw INFO string from VCF.

    Raises:
        ValueError: If any entry is malformed (missing '=' sign).
    """
    for entry in info_str.split(";"):
        if entry and "=" not in entry:
            raise ValueError(
                f"Malformed INFO entry: {entry} (expected 'key=value' pair)"
            )


def parse_format(format_str: str, sample_str: str) -> dict[str, str]:
    """
    Parse the FORMAT and sample fields from a VCF line.

    Args:
        format_str: The FORMAT field from the VCF, e.g. "GT:GQ"
        sample_str: The sample field from the VCF, e.g. "0/1:76"

    Returns:
        Dictionary mapping FORMAT keys to sample values, e.g. {"GT": "0/1", "GQ": "76"}
    """
    return dict(zip(format_str.split(":"), sample_str.split(":")))


def parse_vcf_line(
    line: str,
    vep_fields: list[str],
    format_fields: list[str],
    parse_info: Callable[[str], dict[str, str]],
    parse_format: Callable[[str, str], dict[str, str]],
    validate_info_values: Callable[[dict[str, str]], dict[str, int]],
    validate_info_structure: Callable[[str], None],
    ncol: int = 10,
) -> dict[str, str]:
    """
    Parse a single VCF line into a dictionary.

    Args:
        line: A line from a VCF file.
        vep_fields: List of VEP annotation fields.
        format_fields: List of FORMAT fields to extract.
        parse_info: a function to parse INFO field.
        parse_format: a function to parse FORAMT field.
        validate_info_values: Function to validate INFO parsing.
        validate_info_structure: Function to validate INFO parsing.
        ncol: number of columns from VCF file to process.

    Returns:
        Dictionary with parsed fields.
    """
    columns = line.strip().split("\t")
    if len(columns) < ncol:
        raise ValueError(f"Less than {ncol} columns in VCF line: {line}")
    chrom, pos, _, ref, alt, qual, fltr, info, fmt, sample = columns[:ncol]

    # turn INFO columninto a dict
    # INFO: FAU=46;FCU=28
    # dict: {FAU: 46, FCU: 28}
    info_dict = parse_info(info)
    validate_info_values(info_dict)
    validate_info_structure(info)

    # match values from FORMAT and SAMPLE columns
    format_dict = parse_format(fmt, sample)

    # no default here to easier check if CSQ exists
    csq_data = info_dict.get("CSQ")
    csq_dict = {}
    if csq_data:
        # take the first annotation
        first_annotation = csq_data.split(",")[0]
        csq_values = first_annotation.split("|")
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

    # populate row with values from FORMAT column
    for key in format_fields:
        if key not in format_dict:
            raise KeyError(f"Required FORMAT field '{key}' is missing in sample data.")
        row[key] = format_dict[key]

    return row


def parse_sv_vcf_line(
    line: str,
    parse_info: Callable[[str], dict[str, str]],
    parse_format: Callable[[str, str], dict[str, str]],
    validate_info_values: Callable[[dict[str, str]], dict[str, int]],
    validate_info_structure: Callable[[str], None],
    ncol: int = 10,
) -> dict[str, str]:
    """
    Parse a single structural variant VCF line into a dictionary.

    Args:
        line: A line from a VCF file.
        parse_info: Function to parse INFO field.
        parse_format: Function to parse FORAMT field.
        validate_info_values: Function to validate INFO parsing.
        validate_info_structure: Function to validate INFO parsing.
        ncol: number of columns from the VCF file to process.

    Returns:
        Dictionary with parsed fields
    """
    parts = line.strip().split("\t")
    if len(parts) < ncol:
        raise ValueError(f"Less than {ncol} columns in VCF line: {line}")
    chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample_data = parts[:ncol]

    # Parse ID field
    match = re.search(r"\.(.*?)\.", id_)
    id_clean = match.group(1) if match else id_

    # Parse INFO field
    info_dict = parse_info(info)
    validate_info_structure(info)
    validate_info_values(info_dict)

    # Parse FORMAT and sample data
    format_dict = parse_format(format_, sample_data)

    # this may require subsetting depending on your needs
    row = {
        "CHROM": chrom,
        "POS": pos,
        # available for INS & DEL only
        "END": info_dict.get("END", ""),
        "TYPE": id_clean,
        "SVLEN": info_dict.get("SVLEN", ""),
        "ALT": alt,
        "FILTER": filter_,
        "COVERAGE": info_dict.get("COVERAGE", ""),
        "SUPPORT": info_dict.get("SUPPORT", ""),
        "STRAND": info_dict.get("STRAND", ""),
        "VAF": info_dict.get("VAF", ""),
        "GENOTYPE": format_dict.get("GT", ""),
        "GENOME QUALITY": format_dict.get("GQ", ""),
        "DEPTH REF": format_dict.get("DR", ""),
        "DEPTH TRANS": format_dict.get("DV", ""),
    }
    return row


def safe_convert(value: str, target_type: Callable, default=None):
    """Safely convert a string to a target type, return default on failure."""
    try:
        return target_type(value)
    except (ValueError, TypeError):
        return default


def parse_cnvkit_vcf_line(
    vcf_line: str,
    parse_info: Callable[[str], dict[str, str]],
    parse_format: Callable[[str], dict[str, str]],
    validate_info_values: Callable[[dict[str, str]], dict[str, int]],
    validate_info_structure: Callable[[str], None],
    ncol: int = 10,
) -> dict[str, str]:
    """
    Parse a single CNVkit VCF line into a dictionary.

    Args:
        vcf_line: A line from a CNVkit VCF file.
        parse_info: Function to parse INFO field.
        parse_format: Function to parse FORAMT field.
        validate_info_values: Function to validate INFO parsing.
        validate_info_structure: Function to validate INFO parsing.
        ncol: number of columns in the VCF file

    Returns:
        Dictionary with parsed fields.
    """
    # Split the line into its tab-separated fields
    fields = vcf_line.strip().split("\t")
    if len(fields) < ncol:
        raise ValueError(f"Less than {ncol} columns in VCF line: {vcf_line}")

    # unpack columns
    chrom, pos, _, _, alt, _, _, info, format_str, sample_str = fields[:ncol]

    pos = safe_convert(pos, int, 0)
    variant_type = alt.strip("<>")  # Remove angle brackets from ALT field

    # Parse INFO field
    info_dict = parse_info(info)
    validate_info_structure(info)
    validate_info_values(info_dict)

    # Extract and convert INFO fields
    genes = info_dict.get("Genes", "")
    end = safe_convert(info_dict.get("END", ""), int, 0)
    svlen = safe_convert(info_dict.get("SVLEN", ""), int, 0)
    log_odds_ratio = safe_convert(
        info_dict.get("LOG_ODDS_RATIO", ""), float, float("nan")
    )
    corr_cn = safe_convert(info_dict.get("CORR_CN", ""), float, float("nan"))
    probes = safe_convert(info_dict.get("PROBES", ""), int, 0)
    baf_raw = info_dict.get("BAF", "")
    baf = safe_convert(baf_raw, float, baf_raw)  # fallback to raw if not numeric

    # Parse FORMAT and sample fields
    format_dict = parse_format(format_str, sample_str)
    gt = format_dict.get("GT", "")
    cnq = safe_convert(format_dict.get("CNQ", ""), float, float("nan"))
    dp = safe_convert(format_dict.get("DP", ""), float, float("nan"))

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


@contextmanager
def open_vcf(vcf_path: str) -> TextIO:
    """
    Open a VCF file, handling gzip if necessary.

    Args:
        vcf_path: Path to the VCF file.

    Returns:
        File object in text mode.
    """
    path = Path(vcf_path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def vcf_to_df(
    vcf_path: str,
    vep_fields: list,
    format_fields: list,
    parse_vcf_line: Callable[[str, list, list], dict],
    parse_info: Callable[[str], dict[str, str]],
    parse_format: Callable[[str, str], dict[str, str]],
    validate_info_values: Callable[[dict[str, str]], dict[str, int]],
    validate_info_structure: Callable[[str], None],
    ncol: int = 10,
) -> pd.DataFrame:
    """
    Convert a VEP annotated VCF file to a DataFrame.

    Args:
        vcf_path: Path to the VCF file.
        vep_fields: List of VEP annotation fields.
        format_fields: List of FORMAT fields.
        parse_vcf_line: Function to parse a VCF line.
        parse_info: Function to parse INFO field.
        parse_format: Function to parse FORAMT field.
        validate_info_values: Function to validate INFO parsing.
        validate_info_structure: Function to validate INFO parsing.
        ncol: number of columns in the VCF file

    Returns:
        DataFrame with parsed VCF data
    """
    rows = []
    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            parsed_row = parse_vcf_line(
                line,
                vep_fields,
                format_fields,
                parse_info,
                parse_format,
                validate_info_values,
                validate_info_structure,
                ncol,
            )
            rows.append(parsed_row)
    return pd.DataFrame(rows, columns=vep_fields + format_fields)


def sv_vcf_to_df(
    vcf_path: str,
    parse_sv_vcf_line: Callable[[str], dict],
    parse_cnvkit_vcf_line: Callable[[str], dict],
    parse_info: Callable[[str], dict[str, str]],
    parse_format: Callable[[str], dict[str, str]],
    validate_info_values: Callable[[dict[str, str]], dict[str, int]],
    validate_info_structure: Callable[[str], None],
    ncol: int = 10,
    cnvkit: bool = False,
) -> pd.DataFrame:
    """
    Convert a structural variant VCF file to a DataFrame.

    Args:
        vcf_path: Path to the VCF file.
        parse_sv_vcf_line: Function to parse a VCF line with SVs.
        parse_cnvkit_vcf_line: Function to parse a VCF line with CNVs.
        parse_info: Function to parse INFO field.
        parse_format: Function to parse FORAMT field.
        validate_info_values: Function to validate INFO parsing.
        validate_info_structure: Function to validate INFO parsing.
        ncol: number of columns in the VCF file
        cnvkit: whether to parse CNVkit-specific fields.

    Returns:
        DataFrame with parsed VCF data.
    """
    parse_line = parse_cnvkit_vcf_line if cnvkit else parse_sv_vcf_line

    rows = []
    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            try:
                rows.append(
                    parse_line(
                        line,
                        parse_info,
                        parse_format,
                        validate_info_values,
                        validate_info_structure,
                        ncol,
                    )
                )
            except Exception as e:
                raise ValueError(f"Malformed VCF line: {line}") from e

    return pd.DataFrame(rows)


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
            "Some required parameters are missing. Check your config file!"
        )

    # read SNV vcf file
    logging.info("Reading provided VCF files")
    snv_all_df = vcf_to_df(vcf_snv, vep_fields, format_fields, parse_vcf_line)

    # remove not important SNV categories and those not passing default filter
    snv_all_df = snv_all_df[
        (~snv_all_df["Consequence"].isin(snvs_remove))
        & (snv_all_df["FILTER"] == "PASS")
    ]

    # keep only chosen columns
    snv_picked_columns = snv_all_df[columns_keep]

    # rename SYMBOL to GENE for clarity
    snv_picked_columns = snv_picked_columns.rename(columns={"SYMBOL": "GENE"})

    # Collect TP53 SNV to a separate dataframe
    snv_tp53 = snv_picked_columns[snv_picked_columns["GENE"] == "TP53"]
    logging.info(f"TP53 SNVs after filtering: {len(snv_tp53)}")

    # Collect the rest of SNVs to a separate dataframe
    snv_rest = snv_picked_columns[snv_picked_columns["GENE"] != "TP53"]
    logging.info(f"Not TP53 SNVs after filtering: {len(snv_rest)}")

    # read SV vcf file
    try:
        sv_df = sv_vcf_to_df(
            vcf_sv, parse_sv_vcf_line, parse_cnvkit_vcf_line, cnvkit=False
        )
    except ValueError as e:
        logging.warning(e)

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
    cnv_df = sv_vcf_to_df(
        vcf_cnv, parse_sv_vcf_line, parse_cnvkit_vcf_line, cnvkit=True
    )
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
