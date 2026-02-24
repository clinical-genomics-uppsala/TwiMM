import argparse
import functools
import gzip
import logging
import re
import sys
from pathlib import Path
from typing import Any
from typing import Callable
from typing import TextIO

import pandas as pd
import yaml


# --- Decorators ---


def log_execution(func: Callable) -> Callable:
    """Log the start, end, and errors of a function."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logging.info(f"Started: {func.__name__}")
        try:
            result = func(*args, **kwargs)
            logging.info(f"Finished: {func.__name__}")
            return result
        except Exception as e:
            logging.error(f"Error in {func.__name__}: {e}")
            raise

    return wrapper


def safe_parser(func: Callable) -> Callable:
    """
    Decorator for row parsing functions.
    Wraps the function to catch parsing errors and provide context.
    """

    @functools.wraps(func)
    def wrapper(line: str, *args, **kwargs):
        try:
            return func(line, *args, **kwargs)
        except (ValueError, KeyError) as e:
            raise ValueError(f"Parsing error in line '{line.strip()}': {e}") from e

    return wrapper


# --- Helper Functions ---


def parse_info(info_str: str) -> dict[str, str]:
    """
    Parse the INFO field from a VCF line

    Args:
        info_str: The INFO field from the VCF, e.g. "FAU=46;FCU=28"

    Returns:
        Dictionary mapping INFO keys to values.
    """
    return dict(entry.split("=", 1) for entry in info_str.split(";") if "=" in entry)


def parse_format(format_str: str, sample_str: str) -> dict[str, str]:
    """
    Parse the FORMAT and SAMPLE fields from a VCF line.

    Args:
        format_str: The FORMAT field from the VCF, e.g. "GT:GQ"
        sample_str: The sample field from the VCF, e.g. "0/1:76"

    Returns:
        Dictionary mapping FORMAT keys to sample values.
    """
    format_list = format_str.split(":")
    sample_list = sample_str.split(":")
    if len(format_list) != len(sample_list):
        raise ValueError("Sample and format fields have different length and cannot be matched!")
    return dict(zip(format_list, sample_list))


def safe_convert(value: str, target_type: Callable, default=None):
    """Safely convert a string to a target type, return default on failure."""
    try:
        return target_type(value)
    except (ValueError, TypeError):
        return default


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


def get_genes_from_bed(bed_path: str) -> set[str]:
    """
    Extract gene names from the 4th column of a BED file.
    """
    genes = set()
    if not bed_path or not Path(bed_path).exists():
        logging.warning(f"BED file {bed_path} not found or not specified.")
        return genes

    with open(bed_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                genes.add(parts[3])
    return genes


def write_excel_sheet(writer: pd.ExcelWriter, df: pd.DataFrame, sheet_name: str):
    """
    Write a dataframe to an Excel sheet and apply autofilter to the header.
    """
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    worksheet = writer.sheets[sheet_name]

    # Calculate column widths based on maximum length of content or header
    for i, col in enumerate(df.columns):
        # Find maximum length of data in the column
        column_data = df[col].astype(str).str.len()
        column_len = column_data.max()
        # Handle empty DataFrames or all-NaN columns where max() returns NaN or is not computable
        if pd.isna(column_len):
            column_len = 0
        # Compare with the length of the column header itself
        max_len = max(float(column_len), float(len(str(col)))) + 2
        # Set the column width
        worksheet.set_column(i, i, max_len)

    max_row, max_col = df.shape
    if max_row > 0:
        worksheet.autofilter(0, 0, max_row, max_col - 1)


# --- Parsers ---


@safe_parser
def parse_vcf_line(
    line: str,
    vep_fields: list[str],
    format_fields: list[str],
    ncol: int = 10,
) -> dict[str, str]:
    """
    Parse a single VCF line into a dictionary.
    """
    columns = line.strip().split("\t")
    if len(columns) < ncol:
        raise ValueError(f"Less than {ncol} columns in VCF line: {line}")
    chrom, pos, _, ref, alt, qual, fltr, info, fmt, sample = columns[:ncol]

    # turn INFO column into a dict
    info_dict = parse_info(info)

    # match values from FORMAT and SAMPLE columns
    format_dict = parse_format(fmt, sample)

    # no default here to easier check if CSQ exists
    csq_data = info_dict.get("CSQ")
    csq_dict = {}
    if csq_data:
        # take the first annotation
        first_annotation = csq_data.split(",")[0]
        csq_values = first_annotation.split("|")
        # check if you want to collect more VEP values than actually exist
        if len(vep_fields) > len(csq_values):
            raise ValueError(f"Too may VEP fields are requested: VCF line doesn't have that many")
        csq_dict = dict(zip(vep_fields, csq_values))
    else:
        raise ValueError("Could not parse CSQ field. Does it exist?")

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


@safe_parser
def parse_sv_vcf_line(
    line: str,
    allowed_var_types: list[str] = ["DEL", "INS", "INV", "DUP", "BND"],
    info_values: list[str] = ["END", "SVLEN", "COVERAGE", "SUPPORT", "STRAND", "VAF"],
    ncol: int = 10,
) -> dict[str, str]:
    """
    Parse a single structural variant VCF line into a dictionary.
    """
    parts = line.strip().split("\t")
    if len(parts) < ncol:
        raise ValueError(f"Less than {ncol} columns in VCF line: {line}")
    chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample_ = parts[:ncol]

    # Parse ID field
    match = re.search(r"\.(.*?)\.", id_)
    var_type = match.group(1) if match else id_
    # check that id_ is not empty
    if not id_:
        raise ValueError("ID field (type of variant) is empty!")
    # check if the extracted variant type is valid
    if var_type not in allowed_var_types:
        raise ValueError(f"{var_type} not in allowed variants which are: {allowed_var_types}")

    # Parse INFO field
    info_dict = parse_info(info)
    # check if it is empty
    if not info_dict:
        raise ValueError("Could not parse INFO field. Does it exist?")
    # check if all keys exist
    if var_type == "BND":
        # exclude END & SVLEN
        for k in [value for value in info_values if value != "END" and value != "SVLEN"]:
            if k not in info_dict:
                raise ValueError(f"{k} not found in the INFO field!")
    else:
        for k in info_values:
            if k not in info_dict:
                raise ValueError(f"{k} not found in the INFO field!")

    # Parse FORMAT and SAMPLE columns
    format_dict = parse_format(format_, sample_)
    for k in ["GT", "GQ", "DR", "DV"]:
        if k not in format_dict:
            raise ValueError(f"{k} not found in the FORMAT field!")

    # this may require subsetting depending on your needs
    row = {
        "CHROM": chrom,
        "POS": pos,
        # available for INS & DEL only
        "END": info_dict.get("END", ""),
        "TYPE": var_type,
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


@safe_parser
def parse_cnvkit_vcf_line(
    vcf_line: str,
    ncol: int = 10,
) -> dict[str, str]:
    """
    Parse a single CNVkit VCF line into a dictionary.
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
    if not info_dict:
        raise ValueError("Could not parse INFO field!")
    required_info_keys = ["SVLEN", "END", "LOG_ODDS_RATIO", "CORR_CN", "PROBES", "BAF"]
    if "Genes" not in info_dict:
        logging.warning("INFO field missing 'Genes' for CNV entry at %s:%s; continuing with empty GENE field.", chrom, pos)
        genes = ""
    else:
        genes = info_dict.get("Genes", "")

    # Ensure all other required keys are present
    for k in required_info_keys:
        if k not in info_dict:
            raise ValueError(f"'{k}' not found in the INFO field!")

    # Extract and convert INFO fields
    # Genes have been extracted above
    end = safe_convert(info_dict.get("END", ""), int, 0)
    svlen = safe_convert(info_dict.get("SVLEN", ""), int, 0)
    log_odds_ratio = safe_convert(info_dict.get("LOG_ODDS_RATIO", ""), float, float("nan"))
    corr_cn = safe_convert(info_dict.get("CORR_CN", ""), float, float("nan"))
    probes = safe_convert(info_dict.get("PROBES", ""), int, 0)
    baf_raw = info_dict.get("BAF", "")
    baf = safe_convert(baf_raw, float, baf_raw)  # fallback to raw if not numeric

    # Parse FORMAT and sample fields
    format_dict = parse_format(format_str, sample_str)
    for k in ["GT", "CN", "CNQ", "DP"]:
        if k not in format_dict:
            raise ValueError(f"{k} not found in the FORMAT field!")

    # Extract and convert FORMAT-SAMPLE values
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


# --- Main Processing Logic ---


@log_execution
def process_vcf(vcf_path: str, line_parser: Callable[..., dict], **kwargs) -> pd.DataFrame:
    """
    Generic function to process a VCF file into a DataFrame using a specific line parser.

    Args:
        vcf_path: Path to the VCF file.
        line_parser: Function to parse each line of the VCF.
        **kwargs: Additional arguments to pass to the line_parser.

    Returns:
        DataFrame containing the parsed data.
    """
    rows = []
    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            parsed_row = line_parser(line, **kwargs)
            rows.append(parsed_row)
    return pd.DataFrame(rows)


def create_report(snakemake_obj: Any):
    """
    Main pipeline execution function using snakemake object.
    """
    # Set up logging
    log_file = snakemake_obj.log[0] if snakemake_obj.log else None
    logging.basicConfig(
        filename=log_file,
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
        level=logging.INFO,
        filemode="w",
    )

    logging.info("Script started")
    logging.info(f"Sample name: {snakemake_obj.wildcards.sample}")

    # Get input and output paths
    vcf_snv = snakemake_obj.input.vcf_snv
    vcf_sv = snakemake_obj.input.vcf_sv
    vcf_cnv = snakemake_obj.input.vcf_cnv
    output_xlsx = snakemake_obj.output.xlsx

    logging.info(f"Input files: SNV VCF: {vcf_snv}, SV VCF: {vcf_sv}, CNV VCF: {vcf_cnv}\nOutput file: {output_xlsx}")

    # Get params
    filter_yaml_file = snakemake_obj.params.filter_config
    with open(filter_yaml_file) as file:
        filters = yaml.load(file, Loader=yaml.FullLoader)

    format_fields = filters.get("format_fields", [])
    vep_fields = filters.get("vep_info_fields", [])
    columns_keep_snv = filters.get("columns_keep_snv", [])
    snvs_remove = filters.get("snvs_remove", [])
    columns_drop_tn = filters.get("columns_drop_tn", [])
    # keep legacy name idid for clarity
    # idid stands for Insertion, Deletion, Duplication, Inversion (SV)
    idid_min_len = filters.get("sv_min_len", 1000)
    columns_drop_idid = filters.get("columns_drop_sv", [])

    genes_bed = getattr(snakemake_obj.params, "genes_bed", "")
    genes_to_keep = get_genes_from_bed(genes_bed)
    logging.info(f"Loaded {len(genes_to_keep)} genes from {genes_bed}")

    if any(
        x is None
        for x in [snvs_remove, format_fields, vep_fields, columns_keep_snv, idid_min_len, columns_drop_tn, columns_drop_idid]
    ):
        logging.error("Missing parameters")
        raise ValueError("Some required parameters are missing. Check your config file!")

    # --- 1. Process SNV VCF ---
    logging.info(f"Parsing file: {vcf_snv}")
    try:
        snv_all_df = process_vcf(vcf_snv, parse_vcf_line, vep_fields=vep_fields, format_fields=format_fields)
    except Exception as e:
        logging.error(f"{e}")
        snv_all_df = pd.DataFrame()

    logging.info(f"Total number of rows: {len(snv_all_df)}.")

    if not snv_all_df.empty:
        logging.info("Filtering out unwanted SNVs and SNVs not passing the tool's default filter.")
        snv_all_df = snv_all_df[(~snv_all_df["Consequence"].isin(snvs_remove)) & (snv_all_df["FILTER"] == "PASS")]
        logging.info(f"Rows left: {len(snv_all_df)}.")

        logging.info("Removing unwanted columns.")
        # Ensure columns exist before selecting
        existing_cols = [c for c in columns_keep_snv if c in snv_all_df.columns]

        # Move ALT after POS
        if "POS" in existing_cols and "ALT" in existing_cols:
            existing_cols.remove("ALT")
            pos_idx = existing_cols.index("POS")
            existing_cols.insert(pos_idx + 1, "ALT")

        snv_picked_columns = snv_all_df[existing_cols]
        snv_picked_columns = snv_picked_columns.rename(columns={"SYMBOL": "GENE"})

        # Strip transcript prefix from HGVSc and HGVSp
        for col in ["HGVSc", "HGVSp"]:
            if col in snv_picked_columns.columns:
                snv_picked_columns[col] = snv_picked_columns[col].str.replace(r"^.*:", "", regex=True)

        snv_tp53 = snv_picked_columns[snv_picked_columns.get("GENE") == "TP53"]
        snv_rest = snv_picked_columns[snv_picked_columns.get("GENE") != "TP53"]

        if genes_to_keep:
            logging.info(f"Filtering SNV tab to keep only {len(genes_to_keep)} genes.")
            snv_rest = snv_rest[snv_rest["GENE"].isin(genes_to_keep)]
    else:
        snv_tp53 = pd.DataFrame()
        snv_rest = pd.DataFrame()

    logging.info(f"Number of SNVs in TP53: {len(snv_tp53)}")
    logging.info(f"Number of SNVs in the other genes: {len(snv_rest)}")

    # --- 2. Process SV VCF ---
    logging.info(f"Parsing file: {vcf_sv}")
    try:
        sv_df = process_vcf(vcf_sv, parse_sv_vcf_line)
    except Exception as e:
        logging.error(f"{e}")
        sv_df = pd.DataFrame()

    logging.info(f"Total number of rows: {len(sv_df)}")

    # Specific Filter: Translocations from chr4 and chr14
    if not sv_df.empty:
        tn_chr4 = sv_df[(sv_df["CHROM"] == "chr4") & (sv_df["TYPE"] == "BND")]
        tn_chr14 = sv_df[(sv_df["CHROM"] == "chr14") & (sv_df["TYPE"] == "BND")]

        sv_chr14_pass = sv_df[(sv_df["CHROM"] == "chr14") & (sv_df["FILTER"] == "PASS")]
        # Use assignment to avoid SettingWithCopyWarning
        sv_chr14_pass = sv_chr14_pass.copy()
        sv_chr14_pass["SVLEN"] = pd.to_numeric(sv_chr14_pass["SVLEN"], errors="coerce")
        sv_chr14_idid = sv_chr14_pass[(~sv_chr14_pass["TYPE"].isin(["BND"])) & (sv_chr14_pass["SVLEN"].abs() >= idid_min_len)]

        # Drop columns and merge translocations
        tn_chr4 = tn_chr4.drop(columns=[c for c in columns_drop_tn if c in tn_chr4.columns])
        tn_chr14 = tn_chr14.drop(columns=[c for c in columns_drop_tn if c in tn_chr14.columns])
        sv_chr14_idid = sv_chr14_idid.drop(columns=[c for c in columns_drop_idid if c in sv_chr14_idid.columns])

        translocations_df = pd.concat([tn_chr4, tn_chr14], ignore_index=True)
    else:
        translocations_df = pd.DataFrame()
        sv_chr14_idid = pd.DataFrame()

    logging.info(f"Number of translocations from chr4 and chr14: {len(translocations_df)}")
    logging.info(f"Total IDID variants on chr14: {len(sv_chr14_idid)}")

    # --- 3. Process CNV VCF ---
    logging.info(f"Parsing file: {vcf_cnv}")
    try:
        cnv_df = process_vcf(vcf_cnv, parse_cnvkit_vcf_line)
        logging.info("Formatting POS, END, SVLEN in CNV tab to Megabases.")
        # Rename columns to include (Mb)
        cnv_df = cnv_df.rename(columns={"POS": "POS (Mb)", "END": "END (Mb)", "SVLEN": "SVLEN (Mb)"})
        for col in ["POS (Mb)", "END (Mb)", "SVLEN (Mb)"]:
            if col in cnv_df.columns:
                # Convert to numeric and divide by 1e6
                cnv_df[col] = pd.to_numeric(cnv_df[col], errors="coerce") / 1e6
    except Exception as e:
        logging.error(f"{e}")
        cnv_df = pd.DataFrame()

    logging.info(f"Total CNVs read: {len(cnv_df)}")

    # --- Write to Excel ---
    logging.info("Writing XLSX file.")
    with pd.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        write_excel_sheet(writer, snv_rest, "SNV")
        write_excel_sheet(writer, snv_tp53, "TP53")
        write_excel_sheet(writer, translocations_df, "Translocations")
        write_excel_sheet(writer, sv_chr14_idid, "SV")
        write_excel_sheet(writer, cnv_df, "CNV")

    logging.info("Script finished successfully")


if __name__ == "__main__":
    try:
        # Check if 'snakemake' is defined (provided by Snakemake)
        create_report(snakemake)  # type: ignore
    except NameError:
        # If 'snakemake' is not defined, use argparse
        parser = argparse.ArgumentParser(description="Compile XLSX report from VCF files.")
        parser.add_argument("--vcf-snv", required=True, help="Path to SNV VCF file")
        parser.add_argument("--vcf-sv", required=True, help="Path to SV VCF file")
        parser.add_argument("--vcf-cnv", required=True, help="Path to CNV VCF file")
        parser.add_argument("--output-xlsx", "-o", required=True, help="Path to output XLSX file")
        parser.add_argument("--filter-config", required=True, help="Path to filter config YAML file")
        parser.add_argument("--genes-bed", help="Path to genes BED file", default="")
        parser.add_argument("--sample", help="Sample name", default="unknown")
        parser.add_argument("--log", help="Path to log file", default=None)

        args = parser.parse_args()

        # Create a mock object that behaves like the snakemake object
        class SnakemakeMock:
            def __init__(self, args):
                self.input = argparse.Namespace(vcf_snv=args.vcf_snv, vcf_sv=args.vcf_sv, vcf_cnv=args.vcf_cnv)
                self.output = argparse.Namespace(xlsx=args.output_xlsx)
                self.params = argparse.Namespace(filter_config=args.filter_config, genes_bed=args.genes_bed)
                self.wildcards = argparse.Namespace(sample=args.sample)
                self.log = [args.log] if args.log else []

        create_report(SnakemakeMock(args))
