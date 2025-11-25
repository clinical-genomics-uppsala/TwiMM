import pytest
from workflow.scripts.compile_xlsx_report import parse_info, validate_info_values, validate_info_structure


# --- Tests for parse_info ---
def test_parse_info_basic():
    assert parse_info("FAU=46;FCU=28") == {"FAU": "46", "FCU": "28"}

def test_parse_info_empty():
    assert parse_info("") == {}

def test_parse_info_malformed_entries():
    assert parse_info("FAU;FCU=28") == {"FCU": "28"}

# --- Tests for validate_info_values ---
def test_validate_info_values_valid():
    info = {"FAU": "46", "FCU": "28"}
    assert validate_info_values(info) == {"FAU": 46, "FCU": 28}

@pytest.mark.parametrize("info_dict", [
    {"ANN": "upstream_gene_variant"},
    {"X": "0.1"},
    {"Y": "12a"},
    {"Z": ""},
])
def test_validate_info_values_parametrized_invalid(info_dict):
    with pytest.raises(ValueError):
        validate_info_values(info_dict)


# --- Tests for validate_info_structure ---
def test_validate_info_structure_valid():
    validate_info_structure("FAU=46;FCU=28")  # Should not raise


@pytest.mark.parametrize("info_str", [
    "FAU;FCU=28",                # single malformed entry
    "FAU;XYZ;d;d;FCU=28",        # multiple malformed entries
])
def test_validate_info_structure_malformed(info_str):
    with pytest.raises(ValueError) as excinfo:
        validate_info_structure(info_str)
    assert "Malformed INFO entry" in str(excinfo.value)



