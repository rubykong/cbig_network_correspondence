"""Lazy, versioned, content-verified access to the cbig_network_correspondence_data repo.

Files are downloaded on first use (not at import time) and cached locally via
Pooch. The data repo is pinned to a git tag so results are reproducible across
installs; bump DATA_TAG (and regenerate the registry files with
scripts/make_data_registry.py) whenever a new data release is cut.

The cache location defaults to the platform cache dir, overridable with the
CBIG_NETWORK_CORRESPONDENCE_DATA_DIR environment variable -- useful on a
cluster to pre-populate a shared cache on a login node before running jobs on
worker nodes without internet access (see prefetch_all()).
"""
import gzip
import os
import shutil
from pathlib import Path

import pooch

DATA_REPO = "rubykong/cbig_network_correspondence_data"
DATA_TAG = "v1.1.0"
_BASE_URL = f"https://raw.githubusercontent.com/{DATA_REPO}/{DATA_TAG}/"
_REGISTRY_DIR = Path(__file__).parent / "registry"
_CACHE_ROOT = Path(os.environ.get(
    "CBIG_NETWORK_CORRESPONDENCE_DATA_DIR",
    pooch.os_cache("cbig_network_correspondence"),
))

CORE = pooch.create(path=_CACHE_ROOT / "core", base_url=_BASE_URL, registry=None)
CORE.load_registry(_REGISTRY_DIR / "core.txt")

SPIN = pooch.create(path=_CACHE_ROOT / "spin_rotations", base_url=_BASE_URL, registry=None)
SPIN.load_registry(_REGISTRY_DIR / "spin_rotations.txt")


def fetch_file(relative_path):
    """Download (if not already cached) one core data file, return its local path."""
    return CORE.fetch(relative_path)


def fetch_dir(relative_prefix, name_filter=None):
    """
    Download every core-registry file under relative_prefix whose filename also
    contains name_filter (same substring match used historically to pick a given
    atlas's file(s) out of a shared category folder), and return the local
    directory path containing them.
    """
    prefix = relative_prefix.rstrip("/") + "/"
    matches = [
        key for key in CORE.registry
        if key.startswith(prefix) and (name_filter is None or name_filter in Path(key).name)
    ]
    if not matches:
        raise FileNotFoundError(
            f"No registered data files under {relative_prefix!r} matching {name_filter!r}")
    for key in matches:
        CORE.fetch(key)
    return str(Path(CORE.abspath) / relative_prefix)


def fetch_spin_rotations(data_space):
    """
    Download (if not already cached), reassemble, and decompress the spin-rotation
    permutation model for a given surface space. Returns the local .pkl path.
    """
    prefix = f"spin_rotations/{data_space}/"
    final_path = Path(SPIN.abspath) / data_space / "1000_spin_permutations_state0.pkl"
    if final_path.exists():
        return str(final_path)

    part_keys = sorted(
        (key for key in SPIN.registry if key.startswith(prefix) and ".part" in key),
        key=lambda key: int(key.rsplit(".part", 1)[1]),
    )
    if not part_keys:
        raise FileNotFoundError(f"No spin rotation data registered for space {data_space!r}")

    compressed_path = final_path.with_name(final_path.name + ".gz")
    compressed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(compressed_path, "wb") as f_out:
        for key in part_keys:
            part_path = SPIN.fetch(key)
            with open(part_path, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)

    with gzip.open(compressed_path, "rb") as f_in, open(final_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    compressed_path.unlink()

    return str(final_path)


def prefetch_all(include_spin=False):
    """Warm the local cache for every registered core file (and optionally spin
    rotations too). Useful to run once on a machine with internet access before
    running jobs on a cluster node that has none."""
    for key in CORE.registry:
        CORE.fetch(key)
    if include_spin:
        for data_space in {key.split("/")[1] for key in SPIN.registry}:
            fetch_spin_rotations(data_space)
