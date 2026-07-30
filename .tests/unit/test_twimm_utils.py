import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

TEST_DIR = Path(__file__).parent.resolve()
SCRIPT_DIR = TEST_DIR / "../../workflow/scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from twimm_utils import (  # type: ignore
    expand_cfg_placeholders,
    get_local_vcfs_for_svdb_merge,
    get_svdb_merge_priority,
    get_tc_file,
    get_snv_caller_output,
    get_concat_caller_vcfs,
    get_ubam_input,
)


def wildcards(**kwargs):
    return SimpleNamespace(**kwargs)


# --- Tests for expand_cfg_placeholders ---
@pytest.mark.parametrize(
    "obj, cfg_vars, expected",
    [
        ("{{REF_DATA}}/genome.fasta", {"REF_DATA": "/data/ref"}, "/data/ref/genome.fasta"),
        ("no placeholders here", {"REF_DATA": "/data/ref"}, "no placeholders here"),
        (
            ["{{A}}/x", "{{B}}/y"],
            {"A": "1", "B": "2"},
            ["1/x", "2/y"],
        ),
        (
            {"path": "{{A}}/x", "nested": {"path": "{{B}}/y"}},
            {"A": "1", "B": "2"},
            {"path": "1/x", "nested": {"path": "2/y"}},
        ),
        (42, {}, 42),
    ],
)
def test_expand_cfg_placeholders(obj, cfg_vars, expected):
    assert expand_cfg_placeholders(obj, cfg_vars) == expected


def test_expand_cfg_placeholders_missing_key_raises():
    with pytest.raises(KeyError):
        expand_cfg_placeholders("{{MISSING}}/x", {"REF_DATA": "/data/ref"})


# --- Tests for get_local_vcfs_for_svdb_merge ---
SVDB_CONFIG = {
    "svdb_merge": {
        "tc_method": [
            {
                "name": "purecn",
                "priority": "sv1,sv2",
                "sv_caller": [
                    {"caller": "sv1", "vcf": "sv_calling/sv1/{sample}_{type}.vcf"},
                    {"caller": "sv2", "vcf": "sv_calling/sv2/{sample}_{type}.vcf"},
                ],
            }
        ]
    }
}


def test_get_local_vcfs_for_svdb_merge_no_suffix():
    wc = wildcards(tc_method="purecn", sample="sample1", type="T")
    assert get_local_vcfs_for_svdb_merge(wc, SVDB_CONFIG) == [
        "sv_calling/sv1/sample1_T.vcf",
        "sv_calling/sv2/sample1_T.vcf",
    ]


def test_get_local_vcfs_for_svdb_merge_with_suffix():
    wc = wildcards(tc_method="purecn", sample="sample1", type="T")
    assert get_local_vcfs_for_svdb_merge(wc, SVDB_CONFIG, add_suffix=True) == [
        "sv_calling/sv1/sample1_T.vcf:sv1",
        "sv_calling/sv2/sample1_T.vcf:sv2",
    ]


def test_get_local_vcfs_for_svdb_merge_unknown_tc_method_raises():
    wc = wildcards(tc_method="unknown", sample="sample1", type="T")
    with pytest.raises(KeyError):
        get_local_vcfs_for_svdb_merge(wc, SVDB_CONFIG)


def test_get_local_vcfs_for_svdb_merge_missing_config_section_raises():
    with pytest.raises(KeyError):
        get_local_vcfs_for_svdb_merge(wildcards(tc_method="purecn"), {})


# --- Tests for get_svdb_merge_priority ---
def test_get_svdb_merge_priority():
    wc = wildcards(tc_method="purecn")
    assert get_svdb_merge_priority(wc, SVDB_CONFIG) == "sv1,sv2"


def test_get_svdb_merge_priority_unknown_tc_method_raises():
    wc = wildcards(tc_method="unknown")
    with pytest.raises(KeyError):
        get_svdb_merge_priority(wc, SVDB_CONFIG)


# --- Tests for get_tc_file ---
def test_get_tc_file_pathology():
    wc = wildcards(tc_method="pathology", sample="sample1", type="T")
    assert get_tc_file(wc, {"samples": "config/samples.tsv"}) == "config/samples.tsv"


def test_get_tc_file_non_pathology():
    wc = wildcards(tc_method="purecn", sample="sample1", type="T")
    assert get_tc_file(wc, {}) == "cnv_sv/purecn_purity_file/sample1_T.purity.txt"


# --- Tests for get_snv_caller_output ---
def test_get_snv_caller_output_deepsomatic():
    assert get_snv_caller_output(wildcards(), {"use_deepsomatic": True}) == "snv_indels/snv_concat/{sample}_{type}.vcf.gz"


def test_get_snv_caller_output_clairs_to():
    assert get_snv_caller_output(wildcards(), {"use_deepsomatic": False}) == "snv_indels/clairs_to/{sample}_{type}.vcf.gz"


# --- Tests for get_concat_caller_vcfs ---
def test_get_concat_caller_vcfs():
    wc = wildcards(sample="sample1", type="T")
    assert get_concat_caller_vcfs(wc) == [
        "snv_indels/clairs_to/sample1_T.caller_tagged.vcf.gz",
        "snv_indels/deepsomatic_t_only/sample1_T.caller_tagged.vcf.gz",
    ]


# --- Tests for get_ubam_input ---
def test_get_ubam_input():
    # Mirrors the production Illumina MultiIndex built in common.smk
    # (sample, type, flowcell, lane, barcode); get_ubam_input looks up only
    # by (sample, type), so this exercises the partial-index lookup rather
    # than an exact match on every level.
    units = pd.DataFrame(
        {
            "sample": ["sample1", "sample1", "sample2"],
            "type": ["T", "T", "T"],
            "flowcell": ["fc1", "fc1", "fc1"],
            "lane": ["L1", "L2", "L1"],
            "barcode": ["bc1", "bc1", "bc1"],
            "bam": ["a.bam", "b.bam", "c.bam"],
        }
    ).set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()
    wc = wildcards(sample="sample1", type="T")
    assert get_ubam_input(wc, units) == ["a.bam", "b.bam"]
