import pytest
from workflow.scripts.compile_xlsx_report import (
    parse_info,
    parse_vcf_line,
    parse_format,
    validate_info_structure,
    validate_info_values
)


# --- Tests for parse_info ---
def test_parse_info_basic():
    assert parse_info("FAU=46;FCU=28") == {"FAU": "46", "FCU": "28"}

def test_parse_info_empty():
    assert parse_info("") == {}

def test_parse_info_malformed_entries():
    assert parse_info("FAU;FCU=28") == {"FCU": "28"}

# --- Tests for validate_info_values ---
# def test_validate_info_values_valid():
#     info = {"FAU": "46", "FCU": "28"}
#     assert validate_info_values(info) == {"FAU": 46, "FCU": 28}

# @pytest.mark.parametrize("info_dict", [
#     {"ANN": "upstream_gene_variant"},
#     {"X": "0.1"},
#     {"Y": "12a"},
#     {"Z": ""},
# ])
# def test_validate_info_values_parametrized_invalid(info_dict):
#     with pytest.raises(ValueError):
#         validate_info_values(info_dict)


# --- Tests for validate_info_structure ---
# def test_validate_info_structure_valid():
#     validate_info_structure("FAU=46;FCU=28")  # Should not raise


# @pytest.mark.parametrize("info_str", [
#     "FAU;FCU=28",                # single malformed entry
#     "FAU;XYZ;d;d;FCU=28",        # multiple malformed entries
# ])
# def test_validate_info_structure_malformed(info_str):
#     with pytest.raises(ValueError) as excinfo:
#         validate_info_structure(info_str)
#     assert "Malformed INFO entry" in str(excinfo.value)


# --- Test cases ---
def test_parse_vcf_line_valid():
    line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35"
    vep_fields = ["Gene", "Impact"]
    format_fields = ["GT", "DP"]

    result = parse_vcf_line(
        line=line,
        vep_fields=vep_fields,
        format_fields=format_fields,
        parse_info=parse_info,
        parse_format=parse_format,
        validate_info_values=validate_info_values,
        validate_info_structure=validate_info_structure
    )

    assert result["CHROM"] == "chr1"
    assert result["POS"] == "12345"
    assert result["REF"] == "A"
    assert result["ALT"] == "T"
    assert result["GT"] == "0/1"
    assert result["DP"] == "35"
    assert result["Gene"] == "gene1"
    assert result["Impact"] == "impact1"

def test_parse_vcf_line_missing_columns():
    line = "chr1\t12345\t.\tA\tT"  # Too few columns
    with pytest.raises(ValueError):
        parse_vcf_line(
            line=line,
            vep_fields=[],
            format_fields=[],
            parse_info=parse_info,
            parse_format=parse_format,
            validate_info_values=validate_info_values,
            validate_info_structure=validate_info_structure,

        )

def test_parse_vcf_line_missing_format_field():
    line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28\tGT:DP\t0/1:35"
    format_fields = ["GT", "DP", "GQ"]  # GQ missing
    with pytest.raises(KeyError):
        parse_vcf_line(
            line=line,
            vep_fields=[],
            format_fields=format_fields,
            parse_info=parse_info,
            parse_format=parse_format,
            validate_info_values=validate_info_values,
            validate_info_structure=validate_info_structure,
        )

def test_parse_vcf_line_malformed_info():
    def bad_validate_info_structure(info_str):
        raise ValueError("Malformed INFO entry")
    line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU;FCU=28\tGT:DP\t0/1:35"
    with pytest.raises(ValueError):
        parse_vcf_line(
            line=line,
            vep_fields=[],
            format_fields=["GT", "DP"],
            parse_info=parse_info,
            parse_format=parse_format,
            validate_info_values=validate_info_values,
            validate_info_structure=bad_validate_info_structure,
        )

# def test_parse_vcf_line_invalid_info_values(parse_info, parse_format, validate_info_structure):
#     def bad_validate_info_values(info_dict):
#         raise ValueError("Invalid value")
#     line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=abc;FCU=28\tGT:DP\t0/1:35"
#     with pytest.raises(ValueError):
#         parse_vcf_line(
#             line=line,
#             vep_fields=[],
#             format_fields=["GT", "DP"],
#             parse_info=parse_info,
#             parse_format=parse_format,
#             validate_info_values=bad_validate_info_values,
#             validate_info_structure=validate_info_structure,
#         )

