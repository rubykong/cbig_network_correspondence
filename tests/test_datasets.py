"""Unit tests for datasets.py's lookup/matching logic that don't require
network access. CORE/SPIN are pooch registries built from the checked-in
registry/*.txt files at import time -- no fetch happens until .fetch()/
fetch_file()/fetch_dir() is actually called, so registry-membership and
prefix-matching logic can be tested without downloading anything."""

import importlib

import pytest

from cbig_network_correspondence import datasets


def test_core_registry_contains_known_files():
    assert "atlas_list" in datasets.CORE.registry
    assert "atlas_list_fullname" in datasets.CORE.registry
    assert "atlas_config/BL3" in datasets.CORE.registry


def test_spin_registry_has_part_files_per_surface_space():
    part_keys = [k for k in datasets.SPIN.registry if ".part" in k]
    assert part_keys, "expected at least one spin_rotations .part file in the registry"
    assert any(k.startswith("spin_rotations/fsaverage6/") for k in part_keys)
    assert any(k.startswith("spin_rotations/fs_LR_32k/") for k in part_keys)


def test_fetch_dir_raises_before_any_network_call_when_nothing_matches():
    with pytest.raises(FileNotFoundError):
        datasets.fetch_dir("atlases/FSLMNI2mm", name_filter="NOT_A_REAL_ATLAS_NAME")


def test_fetch_spin_rotations_raises_for_unknown_space():
    with pytest.raises(FileNotFoundError):
        datasets.fetch_spin_rotations("NOT_A_REAL_SPACE")


def test_cache_root_respects_env_var_override(tmp_path, monkeypatch):
    monkeypatch.setenv("CBIG_NETWORK_CORRESPONDENCE_DATA_DIR", str(tmp_path))
    try:
        importlib.reload(datasets)
        assert datasets._CACHE_ROOT == tmp_path
    finally:
        monkeypatch.delenv("CBIG_NETWORK_CORRESPONDENCE_DATA_DIR", raising=False)
        importlib.reload(datasets)
