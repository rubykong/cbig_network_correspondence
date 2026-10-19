"""
Compute overlap_results for a newly promoted atlas against every other existing
atlas, in both directions (overlap_results/<new>/<other>.mat and
overlap_results/<other>/<new>.mat), matching the convention every existing atlas
pair already follows.

Run this as the third step, after stage_atlas.sh and promote_atlas.py, and only
once you've checked the promoted atlas looks correct (right Data_Type, right
number of networks read back, sane overlap values for a known comparison, etc.).

Unlike parse_atlas_config.py/promote_atlas.py, this needs the FULL toolbox
dependency stack (brainspace, nibabel, vtk, ...), since it reuses the toolbox's
own compute_overlap()/read_atlas() for the actual Dice computation -- run it in
the same environment you use for the toolbox itself (e.g. Ruby_py3/NCT_env).
This script always imports from THIS checked-out repo's src/ (prepended to
sys.path below), never from a `pip install`-ed released version that might be
sitting in the same environment -- so it exercises your local, in-development
code (compute_overlap_with_atlases.py, datasets.py, etc.), not whatever was last
published to PyPI.

Deliberately operates on the LOCAL cbig_network_correspondence_data working copy
directly -- no network access at all, not even for shared resources. A
just-promoted atlas isn't tagged/pushed yet, so there is nothing for it to fetch
remotely anyway, and even for the OTHER (already-tagged) atlases, this should
still read whatever is on disk right now rather than quietly fetching from
GitHub, in case local files differ from what's tagged. It reuses the toolbox's
read_atlas()/compute_overlap() functions (the actual Dice math) but builds its
own lightweight LocalAtlasParams instead of RefAtlasParams/AtlasParams, and
monkeypatches datasets.fetch_file()/fetch_dir() (which read_atlas() calls
internally for cortical masks and network-assignment mappings) to resolve
straight from this local checkout instead of downloading -- read_atlas() has no
other way to reach those two resources, so without this patch they'd still hit
the network even though everything else here is local-only.

Usage: compute_overlap.py <config.ini> [--atlases A,B,C] [--force]
"""
import argparse
import sys
from pathlib import Path

from parse_atlas_config import ReleaseConfig

# Force this local checkout's src/ ahead of anything else on sys.path, so this
# script always exercises local, in-development toolbox code -- never a
# `pip install`-ed released version that happens to be sitting in the same
# environment.
TOOLBOX_SRC = str(Path(__file__).resolve().parents[2] / "src")
if TOOLBOX_SRC not in sys.path:
    sys.path.insert(0, TOOLBOX_SRC)

DATA_REPO = Path("/home/rkong/storage/ruby/data/cbig_network_correspondence_data")


class LocalAtlasParams:
    """Drop-in equivalent of RefAtlasParams/AtlasParams that reads straight from
    a local cbig_network_correspondence_data checkout instead of datasets.py's
    tagged-release fetch -- needed because a just-promoted atlas isn't tagged yet."""

    def __init__(self, grab_info, atlas_name, target_space, data_repo):
        atlas_config_path = data_repo / "atlas_config" / atlas_name
        if not atlas_config_path.is_file():
            raise FileNotFoundError(f"No atlas_config for {atlas_name!r}: {atlas_config_path}")
        self.config = grab_info.read_config(str(atlas_config_path))
        self.data_path = str(data_repo / "atlases" / target_space / self.config.category)


def compute_overlap_matrix(coa, grab_info, data_repo, ref_name, other_name):
    """overlap_results/<ref_name>/<other_name>.mat -- ref's networks as rows,
    other's as columns, matching compute_overlap_data_all's existing convention."""
    ref_params = LocalAtlasParams(grab_info, ref_name, _native_space(grab_info, data_repo, ref_name), data_repo)
    other_params = LocalAtlasParams(grab_info, other_name, ref_params.config.space, data_repo)
    return coa.compute_overlap(ref_params, other_params)


def _native_space(grab_info, data_repo, atlas_name):
    config = grab_info.read_config(str(data_repo / "atlas_config" / atlas_name))
    return config.space


def _patch_datasets_to_local(datasets_module, data_repo):
    """
    read_atlas() calls datasets.fetch_file() directly (for cortical masks and,
    when Data_NetworkAssignment is true, the network-assignment mapping) -- the
    only two resources LocalAtlasParams doesn't already resolve itself. Replace
    both functions so every path this script touches resolves from the local
    checkout, with zero network access, instead of datasets.py's normal
    tagged-release fetch.
    """
    def local_fetch_file(relative_path):
        local_path = data_repo / relative_path
        if not local_path.is_file():
            raise FileNotFoundError(
                f"Expected local file not found: {local_path} (no network fetch is "
                "performed by compute_overlap.py -- everything must already be on disk)"
            )
        return str(local_path)

    def local_fetch_dir(relative_prefix, name_filter=None):
        local_dir = data_repo / relative_prefix
        if not local_dir.is_dir():
            raise FileNotFoundError(
                f"Expected local directory not found: {local_dir} (no network fetch is "
                "performed by compute_overlap.py -- everything must already be on disk)"
            )
        return str(local_dir)

    datasets_module.fetch_file = local_fetch_file
    datasets_module.fetch_dir = local_fetch_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="path to the atlas release config.ini")
    parser.add_argument("--atlases", default=None,
                         help="comma-separated atlas names to compute against (default: everything "
                              "in atlas_list except the new atlas itself)")
    parser.add_argument("--force", action="store_true",
                         help="recompute and overwrite overlap_results that already exist")
    args = parser.parse_args()

    try:
        import cbig_network_correspondence.compute_overlap_with_atlases as coa
        import cbig_network_correspondence.datasets as datasets_module
        import cbig_network_correspondence.grab_data_info as grab_info
        import scipy.io
    except ImportError as exc:
        print(f"compute_overlap.py: could not import the toolbox ({exc}). This script needs "
              "the full dependency stack (brainspace, nibabel, vtk, ...) -- run it in the same "
              "environment you use for the toolbox itself, with the package on your path.",
              file=sys.stderr)
        sys.exit(1)

    # sys.path ordering alone doesn't help if some earlier import already cached a
    # different cbig_network_correspondence in sys.modules (Python checks the
    # cache before consulting sys.path) -- so verify explicitly rather than assume.
    if not coa.__file__.startswith(TOOLBOX_SRC):
        print(f"compute_overlap.py: imported cbig_network_correspondence from "
              f"{coa.__file__!r}, not this local checkout ({TOOLBOX_SRC!r}). "
              "A different install is shadowing it (e.g. already imported earlier "
              "in this process, or sys.modules was pre-populated) -- start a fresh "
              "interpreter, or `pip uninstall cbig_network_correspondence` from "
              "this environment if it's pip-installed elsewhere.",
              file=sys.stderr)
        sys.exit(1)

    _patch_datasets_to_local(datasets_module, DATA_REPO)

    cfg = ReleaseConfig(args.config)
    new_name = cfg.name

    if not (DATA_REPO / "atlas_config" / new_name).is_file():
        print(f"compute_overlap.py: {new_name!r} not found in {DATA_REPO}/atlas_config/ -- "
              "run promote_atlas.py first.", file=sys.stderr)
        sys.exit(1)

    if args.atlases:
        other_names = [a.strip() for a in args.atlases.split(",") if a.strip()]
    else:
        all_names = (DATA_REPO / "atlas_list").read_text().splitlines()
        other_names = [a for a in all_names if a != new_name]

    print(f"Computing overlap_results for {new_name!r} against {len(other_names)} atlas(es): "
          f"{', '.join(other_names)}")

    written, skipped, failed = [], [], []
    for other_name in other_names:
        for ref_name, oth_name in ((new_name, other_name), (other_name, new_name)):
            out_path = DATA_REPO / "overlap_results" / ref_name / f"{oth_name}.mat"
            if out_path.exists() and not args.force:
                print(f"  SKIP {ref_name}/{oth_name}.mat (already exists, use --force to recompute)")
                skipped.append((ref_name, oth_name))
                continue
            try:
                overlap_val = compute_overlap_matrix(coa, grab_info, DATA_REPO, ref_name, oth_name)
            except Exception as exc:  # noqa: BLE001
                print(f"  FAIL {ref_name}/{oth_name}.mat: {exc!r}")
                failed.append((ref_name, oth_name, exc))
                continue
            out_path.parent.mkdir(parents=True, exist_ok=True)
            scipy.io.savemat(str(out_path), {"overlap": overlap_val})
            print(f"  wrote {out_path} (shape={overlap_val.shape})")
            written.append((ref_name, oth_name))

    print()
    print(f"Done: {len(written)} written, {len(skipped)} skipped (already existed), {len(failed)} failed.")
    if failed:
        for ref_name, oth_name, exc in failed:
            print(f"  FAILED {ref_name}/{oth_name}.mat: {exc!r}")
        sys.exit(1)
    print(f"This script did NOT run git. Review {DATA_REPO}/overlap_results/, then git add/commit/push yourself.")


if __name__ == "__main__":
    main()
