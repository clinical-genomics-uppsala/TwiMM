import io
import logging
import textwrap

import pytest
import pysam

import sys
from pathlib import Path

TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from fix_af import modifyHeader, writeNewVcf  # type: ignore


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_VCF_HEADER_CLAIRS_TO = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
    """)

_VCF_HEADER_DEEPSOMATIC = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fraction">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
    """)

_VCF_HEADER_NO_AF_VAF = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
    """)


def _write_vcf(tmp_path, name, header_str, data_lines):
    """Write a plain-text VCF and bgzip+tabix it so pysam can read it."""
    vcf_path = tmp_path / name
    vcf_path.write_text(header_str + "".join(data_lines))
    gz_path = tmp_path / (name + ".gz")
    pysam.tabix_compress(str(vcf_path), str(gz_path), force=True)
    pysam.tabix_index(str(gz_path), preset="vcf", force=True)
    return str(gz_path)


def _read_first_record(path):
    vcf = pysam.VariantFile(path)
    records = list(vcf.fetch())
    assert records, "Expected at least one record"
    return records[0]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSnvConcatAfVafFallback:
    """writeNewVcf with caller='snv_concat' AF→VAF fallback logic."""

    def test_clairs_to_record_af_propagated_to_info(self, tmp_path):
        """ClairS-TO style: FORMAT/AF present → INFO/AF set from FORMAT/AF."""
        data = ["chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AF:DP\t0/1:0.35:50\n"]
        in_path = _write_vcf(tmp_path, "clairs_to.vcf", _VCF_HEADER_CLAIRS_TO, data)
        out_path = str(tmp_path / "out_clairs.vcf.gz")

        vcf = pysam.VariantFile(in_path)
        header = modifyHeader("snv_concat", vcf.header)
        writeNewVcf(out_path, header=header, vcf=vcf, caller="snv_concat")

        rec = _read_first_record(out_path)
        assert "AF" in rec.info, "INFO/AF should be set for a ClairS-TO record"
        assert pytest.approx(rec.info["AF"][0]) == 0.35

    def test_deepsomatic_record_vaf_fallback_to_info_af(self, tmp_path):
        """DeepSomatic style: FORMAT/AF absent, FORMAT/VAF present → INFO/AF set from VAF."""
        data = ["chr1\t200\t.\tG\tC\t.\tPASS\t.\tGT:VAF:DP\t0/1:0.42:60\n"]
        in_path = _write_vcf(tmp_path, "deepsomatic.vcf", _VCF_HEADER_DEEPSOMATIC, data)
        out_path = str(tmp_path / "out_deepsomatic.vcf.gz")

        vcf = pysam.VariantFile(in_path)
        header = modifyHeader("snv_concat", vcf.header)
        writeNewVcf(out_path, header=header, vcf=vcf, caller="snv_concat")

        rec = _read_first_record(out_path)
        assert "AF" in rec.info, "INFO/AF should be populated via VAF fallback"
        assert pytest.approx(rec.info["AF"][0]) == 0.42

    def test_neither_af_nor_vaf_emits_warning_and_no_info_af(self, tmp_path, caplog):
        """Record with neither FORMAT/AF nor FORMAT/VAF: INFO/AF absent, warning logged."""
        data = ["chr1\t300\t.\tC\tT\t.\tPASS\t.\tGT:DP\t0/1:30\n"]
        in_path = _write_vcf(tmp_path, "no_af.vcf", _VCF_HEADER_NO_AF_VAF, data)
        out_path = str(tmp_path / "out_no_af.vcf.gz")

        vcf = pysam.VariantFile(in_path)
        header = modifyHeader("snv_concat", vcf.header)

        with caplog.at_level(logging.WARNING, logger="root"):
            writeNewVcf(out_path, header=header, vcf=vcf, caller="snv_concat")

        assert any(
            "neither FORMAT/AF nor FORMAT/VAF" in r.message for r in caplog.records
        ), "Expected a warning about missing AF/VAF"

        rec = _read_first_record(out_path)
        assert "AF" not in rec.info, "INFO/AF must not be set when both AF and VAF are absent"

    def test_no_exception_raised_for_deepsomatic_record(self, tmp_path):
        """writeNewVcf must not raise for a DeepSomatic-style record."""
        data = ["chr1\t400\t.\tT\tA\t.\tPASS\t.\tGT:VAF:DP\t1/1:0.95:80\n"]
        in_path = _write_vcf(tmp_path, "deepsomatic2.vcf", _VCF_HEADER_DEEPSOMATIC, data)
        out_path = str(tmp_path / "out_deepsomatic2.vcf.gz")

        vcf = pysam.VariantFile(in_path)
        header = modifyHeader("snv_concat", vcf.header)
        writeNewVcf(out_path, header=header, vcf=vcf, caller="snv_concat")
