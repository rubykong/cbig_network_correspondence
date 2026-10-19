"""Unit tests for grab_data_info.py -- pure config parsing, no I/O beyond a
local config file, no network. Covers the string/threshold parsing helpers and
read_config()'s defaulting behavior, since several release-pipeline scripts
(scripts/release_atlas/) depend on these exact semantics."""

import math

import pytest

from cbig_network_correspondence import grab_data_info as gdi


@pytest.mark.parametrize(
    "value,expected",
    [
        (None, False),
        ("yes", True),
        ("YES", True),
        ("true", True),
        ("True", True),
        ("t", True),
        ("1", True),
        ("no", False),
        ("false", False),
        ("0", False),
        ("", False),
    ],
)
def test_str2bool(value, expected):
    assert gdi.str2bool(value) is expected


def test_str2num_none():
    assert gdi.str2num(None) is None


def test_str2num_inf():
    assert gdi.str2num("Inf") == math.inf
    assert gdi.str2num("inf") == math.inf
    assert gdi.str2num("-Inf") == -math.inf


def test_str2num_int_and_float():
    assert gdi.str2num("3") == 3
    assert isinstance(gdi.str2num("3"), int)
    assert gdi.str2num("1.5") == 1.5
    assert isinstance(gdi.str2num("1.5"), float)


def test_thresh_str_to_num_none():
    assert gdi.thresh_str_to_num(None) is None


def test_thresh_str_to_num_parses_bracketed_pair():
    assert gdi.thresh_str_to_num("[1.5,Inf]") == (1.5, math.inf)


def test_thresh_str_to_num_tolerates_space_after_comma():
    # the release config format documents "[0, Inf]" (space after comma) as
    # valid -- confirms replace(',', ' ') + split() handles it without
    # producing an extra empty token
    assert gdi.thresh_str_to_num("[0, Inf]") == (0, math.inf)


def _write_config(tmp_path, body):
    config_path = tmp_path / "config.ini"
    config_path.write_text(body)
    return str(config_path)


def test_read_config_hard_atlas_minimal(tmp_path):
    config_path = _write_config(
        tmp_path,
        "[data_info]\n"
        "Data_Category: TestCat\n"
        "Data_Name: TESTATLAS\n"
        "Data_Space: FSLMNI2mm\n"
        "Data_Type: Hard\n",
    )
    args = gdi.read_config(config_path)
    assert args.name == "TESTATLAS"
    assert args.space == "FSLMNI2mm"
    assert args.type == "Hard"
    # netassign wasn't given -- defaults to False
    assert args.netassign is False
    # threshold defaulting only applies to Soft -- Hard stays None
    assert args.threshold is None


def test_read_config_soft_atlas_defaults_threshold(tmp_path):
    config_path = _write_config(
        tmp_path,
        "[data_info]\n"
        "Data_Category: TestCat\n"
        "Data_Name: TESTSOFT\n"
        "Data_Space: FSLMNI2mm\n"
        "Data_Type: Soft\n",
    )
    args = gdi.read_config(config_path)
    assert args.threshold == (0, math.inf)


def test_read_config_soft_atlas_explicit_threshold(tmp_path):
    config_path = _write_config(
        tmp_path,
        "[data_info]\n"
        "Data_Category: TestCat\n"
        "Data_Name: TESTSOFT\n"
        "Data_Space: FSLMNI2mm\n"
        "Data_Type: Soft\n"
        "Data_Threshold: [1.5,Inf]\n",
    )
    args = gdi.read_config(config_path)
    assert args.threshold == (1.5, math.inf)


def test_read_config_metric_atlas_no_default_threshold(tmp_path):
    # matches TW_TASK_NETS in cbig_network_correspondence_data: Metric type
    # with no Data_Threshold stays None -- the [0, Inf] default only applies
    # to Soft, not Metric, despite check_args_format's print message wording
    config_path = _write_config(
        tmp_path,
        "[data_info]\n"
        "Data_Category: TestCat\n"
        "Data_Name: TESTMETRIC\n"
        "Data_Space: FSLMNI2mm\n"
        "Data_Type: Metric\n",
    )
    args = gdi.read_config(config_path)
    assert args.threshold is None


def test_read_config_missing_name_raises(tmp_path):
    config_path = _write_config(
        tmp_path,
        "[data_info]\n"
        "Data_Category: TestCat\n"
        "Data_Space: FSLMNI2mm\n"
        "Data_Type: Hard\n",
    )
    with pytest.raises(ValueError, match="Data_Name"):
        gdi.read_config(config_path)


def test_read_config_missing_type_raises(tmp_path):
    config_path = _write_config(
        tmp_path,
        "[data_info]\n"
        "Data_Category: TestCat\n"
        "Data_Name: TESTATLAS\n"
        "Data_Space: FSLMNI2mm\n",
    )
    with pytest.raises(ValueError, match="Data_Type"):
        gdi.read_config(config_path)
