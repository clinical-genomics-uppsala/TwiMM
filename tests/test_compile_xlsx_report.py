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
    ],
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
        ("GT:GQ:BD", "", pytest.raises(ValueError)),
    ],
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


@pytest.mark.parametrize(
    "line, vep_fields, format_fields, expected",
    [
        # Too few columns
        ("chr1\t12345\t.\tA\tT", [], [], pytest.raises(ValueError)),
        # GQ missing in the line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP", "GQ"],
            pytest.raises(KeyError),
        ),
        # Substitution is missing in the VCF line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact", "Substitution"],
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
        # line without CSQ field
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
    ],
)
def test_parse_vcf_line(line, vep_fields, format_fields, expected):
    with expected:
        parse_vcf_line(line, vep_fields, format_fields, parse_info, parse_format)

# TODO
# --- Tests for parse_sv_vcf_line ---
