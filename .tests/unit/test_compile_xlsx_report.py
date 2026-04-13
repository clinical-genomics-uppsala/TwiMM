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
