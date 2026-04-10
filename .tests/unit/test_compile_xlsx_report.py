import pytest # type: ignore
from pathlib import Path
import sys
import pandas as pd


TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
RESULTS_DIR = TEST_DIR / "../../results"
sys.path.insert(0, str(SCRIPT_DIR))

from compile_xlsx_report import ( # type: ignore
    parse_info,
    parse_vcf_line,
    parse_format,
    parse_cnvkit_vcf_line,
    parse_version_from_container,
    get_genes_from_bed,
    process_sv_vcf,
    filter_cnv_df,
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
        # Short input strings (one item in each)
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
def test_parse_format(format_str, sample_str, expected):
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
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:25",
            ["Gene", "Impact"],
            ["GT", "DP", "GQ"],
            pytest.raises(ValueError),
        ),
        # You request too many VEP fields from the VCF file
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;CSQ=gene1|impact1\tGT:DP\t0/1:35",
            ["Gene", "Impact", "Substitution"],  # Substitution is missing
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
        # CSQ field is missing in VCF line
        (
            "chr1\t12345\t.\tA\tT\t.\tPASS\tFAU=46;FCU=28;\tGT:DP\t0/1:45",
            ["Gene", "Impact"],
            ["GT", "DP"],
            pytest.raises(ValueError),
        ),
    ],
)
def test_parse_vcf_line(line, vep_fields, format_fields, expected):
    if isinstance(expected, dict):
        assert (
            parse_vcf_line(line, vep_fields, format_fields)
            == expected
        )
    else:
        with expected:
            parse_vcf_line(line, vep_fields, format_fields)



@pytest.mark.parametrize(
    "line, expected",
    [
        # Normal case
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=gene1,gene2;SVTYPE=COPY_NORMAL;END=125;SVLEN=124;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            {
                "CHROM": "chr1",
                "POS": 10,
                "VARIANT_TYPE": "COPY_NORMAL",
                "GENE": "gene1,gene2",
                "END": 125,
                "SVLEN": 124,
                "LOG_ODDS_RATIO": -0.1,
                "CORR_CN": 2.0,
                "PROBES": 2,
                "BAF": 0.3,
                "GT": "0/0",
                "CNQ": 80,
                "DP": 3.5,
            },
        ),
        # Genes in INFO is missing (it's OK when the variant is DEL)
        (
            "chr1\t10\t.\tN\t<DEL>\t.\t.\tSVTYPE=COPY_NORMAL;END=125;SVLEN=124;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            {
                "CHROM": "chr1",
                "POS": 10,
                "VARIANT_TYPE": "DEL",
                # Genes missing → empty string
                "GENE": "",
                "END": 125,
                "SVLEN": 124,
                "LOG_ODDS_RATIO": -0.1,
                "CORR_CN": 2.0,
                "PROBES": 2,
                "BAF": 0.3,
                "GT": "0/0",
                "CNQ": 80,
                "DP": 3.5,
            },
        ),
        # Line is too short
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t0/0:2.0:80:3.5:0.3",
            pytest.raises(ValueError),
        ),
        # NO info field
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\t?\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            pytest.raises(ValueError),
        ),
        # One of the expected values in INFO is missing (SVLEN, SVTYPE etc)
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=gene1,gene2;SVTYPE=COPY_NORMAL;END=125;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CN:CNQ:DP:BAF\t0/0:2.0:80:3.5:0.3",
            pytest.raises(ValueError),
        ),
        # One of the expected values in FORMAT is missing (GT, GQ etc)
        (
            "chr1\t10\t.\tN\t<COPY_NORMAL>\t.\t.\tGenes=gene1,gene2;SVTYPE=COPY_NORMAL;END=125;SVLEN=124;LOG_ODDS_RATIO=-0.1;CALLER=cnvkit;CN=NA;CORR_CN=2.0;PROBES=2;BAF=0.3\tGT:CNQ\t0/0:2.0",
            pytest.raises(ValueError),
        ),
    ],
)
def test_parse_cnvkit_vcf_line(line, expected):
    if isinstance(expected, dict):
        assert parse_cnvkit_vcf_line(line) == expected
    else:
        with expected:
            parse_cnvkit_vcf_line(line)


# --- Tests for parse_version_from_container ---
@pytest.mark.parametrize(
    "container, expected",
    [
        ("docker://hydragenetics/severus:1.5", "1.5"),
        ("docker://hydragenetics/pbmm2:1.16", "1.16"),
        ("docker://hkubal/clairs-to:latest", "latest"),
        ("docker://org/tool:", ""),
        ("someimage", "unknown"),
        ("", "unknown"),
    ],
)
def test_parse_version_from_container(container, expected):
    assert parse_version_from_container(container) == expected


# --- Integration test for process_sv_vcf ---
COLO829_VCF = RESULTS_DIR / "COLO829_T.vcf"


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_columns():
    """process_sv_vcf must return all expected columns including the four population frequency columns."""
    df = process_sv_vcf(str(COLO829_VCF))
    expected_cols = [
        "CHROM", "POS", "END", "TYPE", "SVLEN", "ALT", "FILTER", "CALLER",
        "COVERAGE", "SUPPORT", "STRAND", "VAF", "GENOTYPE", "GENOME QUALITY",
        "DEPTH REF", "DEPTH TRANS",
        "GNOMAD_AC", "GNOMAD_AF", "CUSTOM_AC", "CUSTOM_AF",
    ]
    assert list(df.columns) == expected_cols
    assert len(df) > 0


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_types():
    """Population frequency columns must be numeric (not string) where present.
    pandas coerces int+None columns to float64, which is fine for Excel rendering."""
    df = process_sv_vcf(str(COLO829_VCF))
    for col in ["GNOMAD_AC", "GNOMAD_AF", "CUSTOM_AC", "CUSTOM_AF"]:
        assert pd.api.types.is_numeric_dtype(df[col]), \
            f"Column {col} should be numeric, got dtype {df[col].dtype}"
    # Spot-check: at least one row has a non-null GNOMAD_AC value
    assert df["GNOMAD_AC"].notna().any(), "Expected at least one non-null GNOMAD_AC"


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_sample_selection():
    """
    Severus records should have depth from the haplotagged sample (DR/DV present).
    PBSV records should have depth from the AD fallback.
    """
    df = process_sv_vcf(str(COLO829_VCF))
    severus_rows = df[df["CALLER"] == "severus"]
    pbsv_rows = df[df["CALLER"] == "pbsv"]
    assert len(severus_rows) > 0, "Expected at least one severus record in fixture"
    assert len(pbsv_rows) > 0, "Expected at least one pbsv record in fixture"
    # Severus: DR and DV are explicit FORMAT fields → non-empty
    assert (severus_rows["DEPTH TRANS"] != "").all(), "Severus DEPTH TRANS should be non-empty"
    # PBSV: depth comes from AD fallback → non-empty
    assert (pbsv_rows["DEPTH TRANS"] != "").all(), "PBSV DEPTH TRANS should be non-empty (AD fallback)"


@pytest.mark.skipif(not COLO829_VCF.exists(), reason="COLO829_T.vcf not present")
def test_process_sv_vcf_severus_strand_vaf():
    """Severus records should have STRANDS → STRAND populated and VAF from FORMAT."""
    df = process_sv_vcf(str(COLO829_VCF))
    severus_rows = df[df["CALLER"] == "severus"]
    assert len(severus_rows) > 0, "Expected at least one severus record in fixture"
    assert (severus_rows["STRAND"] != "").all(), "Severus STRAND should be non-empty (from STRANDS field)"
    assert (severus_rows["VAF"] != "").all(), "Severus VAF should be non-empty (from FORMAT)"


# --- Tests for get_genes_from_bed ---
def test_get_genes_from_bed(tmp_path):
    # Test valid BED file
    bed_content = (
        "chr1\t100\t200\tGENE1\n"
        "chr2\t300\t400\tGENE2\tignore\n"
        "#comment line\n"
        "\n"
        "chr3\t500\t600\tGENE3\n"
    )
    bed_file = tmp_path / "test.bed"
    bed_file.write_text(bed_content)
    
    genes = get_genes_from_bed(str(bed_file))
    assert genes == {"GENE1", "GENE2", "GENE3"}

    # Test missing file
    assert get_genes_from_bed("non_existent.bed") == set()

    # Test empty path
    assert get_genes_from_bed("") == set()

    # Test malformed line (fewer than 4 columns)
    malformed_content = "chr1\t100\t200\n"
    malformed_file = tmp_path / "malformed.bed"
    malformed_file.write_text(malformed_content)
    assert get_genes_from_bed(str(malformed_file)) == set()


# --- Tests for filter_cnv_df ---

def _make_cnv_df(rows):
    """Build a minimal CNV DataFrame as produced after Mb conversion."""
    return pd.DataFrame(rows, columns=["VARIANT_TYPE", "CORR_CN", "BAF", "SVLEN (Mb)"])


DEFAULTS = dict(corr_cn_min=1.4, corr_cn_max=2.5, baf_min=0.3, baf_max=0.7, svlen_min_mb=0.0)


def test_filter_cnv_df_empty():
    """Empty DataFrame passes through without error."""
    df = pd.DataFrame(columns=["VARIANT_TYPE", "CORR_CN", "BAF", "SVLEN (Mb)"])
    result = filter_cnv_df(df, **DEFAULTS)
    assert result.empty


def test_filter_cnv_df_removes_copy_normal():
    df = _make_cnv_df([
        ("COPY_NORMAL", 2.0, 0.5, 1.0),
        ("DEL", 1.0, 0.1, 2.0),
    ])
    result = filter_cnv_df(df, **DEFAULTS)
    assert "COPY_NORMAL" not in result["VARIANT_TYPE"].values
    assert len(result) == 1


def test_filter_cnv_df_corr_cn_inside_range_removed():
    """Rows with CORR_CN within [1.4, 2.5] are removed."""
    df = _make_cnv_df([
        ("DEL", 1.4, 0.1, 1.0),   # boundary — inside → removed
        ("DEL", 2.5, 0.1, 1.0),   # boundary — inside → removed
        ("DEL", 2.0, 0.1, 1.0),   # inside → removed
        ("DEL", 1.39, 0.1, 1.0),  # outside → kept
        ("DEL", 2.51, 0.1, 1.0),  # outside → kept
    ])
    result = filter_cnv_df(df, **DEFAULTS)
    assert len(result) == 2
    assert list(result["CORR_CN"]) == [1.39, 2.51]


def test_filter_cnv_df_baf_inside_range_removed():
    """Rows with BAF within [0.3, 0.7] are removed."""
    df = _make_cnv_df([
        ("DEL", 0.5, 0.3, 1.0),   # BAF boundary → removed
        ("DEL", 0.5, 0.7, 1.0),   # BAF boundary → removed
        ("DEL", 0.5, 0.5, 1.0),   # BAF inside → removed
        ("DEL", 0.5, 0.29, 1.0),  # BAF outside → kept
        ("DEL", 0.5, 0.71, 1.0),  # BAF outside → kept
    ])
    # CORR_CN=0.5 is outside [1.4, 2.5] so CORR_CN filter passes all rows
    result = filter_cnv_df(df, **DEFAULTS)
    assert len(result) == 2
    assert list(result["BAF"]) == [0.29, 0.71]


def test_filter_cnv_df_nan_corr_cn_removed():
    """Rows with non-numeric CORR_CN are dropped (NaN fails both comparisons)."""
    df = _make_cnv_df([
        ("DEL", ".", 0.1, 1.0),
        ("DEL", None, 0.1, 1.0),
        ("DEL", 1.0, 0.1, 1.0),  # valid, kept
    ])
    result = filter_cnv_df(df, **DEFAULTS)
    assert len(result) == 1
    assert result.iloc[0]["CORR_CN"] == 1.0


def test_filter_cnv_df_nan_baf_removed():
    """Rows with non-numeric BAF are dropped."""
    df = _make_cnv_df([
        ("DEL", 1.0, "NA", 1.0),
        ("DEL", 1.0, 0.1, 1.0),  # valid, kept
    ])
    result = filter_cnv_df(df, **DEFAULTS)
    assert len(result) == 1


def test_filter_cnv_df_svlen_min_zero_skips_filter():
    """svlen_min_mb=0 means no SVLEN filtering is applied."""
    df = _make_cnv_df([
        ("DEL", 1.0, 0.1, 0.0),
        ("DEL", 1.0, 0.1, 0.001),
    ])
    result = filter_cnv_df(df, **{**DEFAULTS, "svlen_min_mb": 0.0})
    assert len(result) == 2


def test_filter_cnv_df_svlen_min_applied():
    """svlen_min_mb > 0 removes rows where |SVLEN (Mb)| < threshold."""
    df = _make_cnv_df([
        ("DEL", 1.0, 0.1, -2.0),  # abs = 2.0 >= 1.0 → kept
        ("DEL", 1.0, 0.1, 0.5),   # abs = 0.5 < 1.0 → removed
        ("DEL", 1.0, 0.1, 1.0),   # abs = 1.0 >= 1.0 → kept
    ])
    result = filter_cnv_df(df, **{**DEFAULTS, "svlen_min_mb": 1.0})
    assert len(result) == 2
    assert list(result["SVLEN (Mb)"]) == [-2.0, 1.0]
