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
@pytest.mark.parametrize(
    "line, vep_fields, format_fields, expected",
    [
        # Normal case
        (
            "chr1\t12345\t.\tA\tT\t33\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP"],
            {
                "CHROM": "chr1",
                "POS": "12345",
                "REF": "A",
                "ALT": "T",
                "QUAL": "33",
                "FILTER": "PASS",
                "GT": "0/1",
                "DP": "35",
                "Gene": "gene1",
                "Impact": "impact1",
            },
        ),
        # Too few columns in VCF line
        ("chr1\t12345\t.\tA\tT", [], [], pytest.raises(ValueError)),
        # GQ missing in the line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP", "GQ"],
            pytest.raises(KeyError),
        ),
        # Substitution is missing in VCF line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact", "Substitution"],
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
        # CSQ field is missing in VCF line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;\tGT:DP\t0/1:35",
            ["Gene", "Impact"],
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
    ],
)
def test_parse_vcf_line(line, vep_fields, format_fields, expected):
    if isinstance(expected, dict):
        assert parse_vcf_line(line, vep_fields, format_fields, parse_info, parse_format) == expected
    else:
        with expected:
            parse_vcf_line(line, vep_fields, format_fields, parse_info, parse_format)
# TODO
# --- Tests for parse_sv_vcf_line ---
