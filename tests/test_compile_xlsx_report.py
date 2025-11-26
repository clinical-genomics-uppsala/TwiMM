import pytest
from workflow.scripts.compile_xlsx_report import (
    parse_info,
    parse_vcf_line,
    parse_format,
)


# --- Tests for parse_info ---
@pytest.mark.parametrize(
        "info_str, expected",
        [
            # Normal case
            ("FAU=46;FCU=28", {"FAU": "46", "FCU": "28"}),

            # empty input
            ("", {}),

            # malformed input
            ("FAU;FCU=28", {"FCU": "28"}),
        ]
)

def test_parse_info(info_str, expected):
    assert parse_info(info_str) == expected


# --- Tests for parse_format ---
@pytest.mark.parametrize(
    "format_str, sample_str, expected",
    [
        # Normal case
        ("GT:GQ", "0/1:76", {"GT": "0/1", "GQ": "76"}),

        # Short input strins (one item in each)
        ("GT", "0/1", {"GT": "0/1"}),

        # Both empty strings
        ("", "", {"": ""}),

        # Non-string inputs → should raise AttributeError because .split() fails
        (None, "0/1:76", pytest.raises(AttributeError)),
        ("GT:GQ", None, pytest.raises(AttributeError)),
        (123, "0/1:76", pytest.raises(AttributeError)),
        ("GT:GQ", 456, pytest.raises(AttributeError)),
        (None, None, pytest.raises(AttributeError)),
        (123, 123, pytest.raises(AttributeError)),
        
        # Unequal length of input strings
        ("GT:GQ:BD", "0/1:76", pytest.raises(ValueError)),
        # one input is empty the other is not
        ("", "0/1:76", pytest.raises(ValueError)),
        ("GT:GQ:BD", "", pytest.raises(ValueError))
    ]
)

def test_parse_format_edge_cases(format_str, sample_str, expected):
    if isinstance(expected, dict):
        assert parse_format(format_str, sample_str) == expected
    else:
        with expected:
            parse_format(format_str, sample_str)



# --- Tests for parse_vcf_line ---
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
            ncol=10,
        )

def test_parse_vcf_line_missing_format_field():
    line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35"
    format_fields = ["GT", "DP", "GQ"]  # GQ missing in the line
    vep_fields = ["Gene", "Impact"]
    with pytest.raises(KeyError):
        parse_vcf_line(
            line=line,
            vep_fields=vep_fields,
            format_fields=format_fields,
            parse_info=parse_info,
            parse_format=parse_format,
        )

def test_parse_vcf_line_missing_vep_field():
    line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35"
    format_fields = ["GT", "DP"]
    vep_fields = ["Gene", "Impact", "Substitution"] # Substitution is missing in the line
    with pytest.raises(ValueError):
        parse_vcf_line(
            line=line,
            vep_fields=vep_fields,
            format_fields=format_fields,
            parse_info=parse_info,
            parse_format=parse_format,
        )

def test_parse_vcf_line_no_csq_field():
    # line without CSQ field
    line = "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;\tGT:DP\t0/1:35"
    vep_fields=["Gene", "Impact"]
    format_fields=["GT", "DP"]
    with pytest.raises(ValueError):
        parse_vcf_line(
            line=line,
            vep_fields=vep_fields,
            format_fields=format_fields,
            parse_info=parse_info,
            parse_format=parse_format,
        )
