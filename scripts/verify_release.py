"""
End-to-end verification to run before cutting a PyPI release.

Exercises the real package (not a reimplementation) against the live, tagged
cbig_network_correspondence_data: one atlas per surface/volume space, the
documented Quick Example, an atlas-vs-atlas comparison, and (optionally) the
spin-rotation permutation path.

Usage (run inside the environment you intend to release with, e.g. Ruby_py3 /
NCT_env -- the full dependency stack must already be installed there, but the
cbig_network_correspondence package itself does not need to be, since this
script always imports it from this checked-out repo's src/, never from a `pip
install`-ed released version that might be sitting in the same environment):

    python scripts/verify_release.py
    python scripts/verify_release.py --skip-spin   # skip the ~1.4GB spin-rotation download/test

See scripts/README.md for what each check means and what to do about failures.
"""
import argparse
import shutil
import sys
import tempfile
import time
import traceback
from pathlib import Path

# Force this local checkout's src/ ahead of anything else on sys.path, so this
# always verifies local, in-development code -- never a `pip install`-ed released
# version that happens to be sitting in the same environment.
TOOLBOX_SRC = str(Path(__file__).resolve().parent.parent / "src")
if TOOLBOX_SRC not in sys.path:
    sys.path.insert(0, TOOLBOX_SRC)

# one representative atlas per space, chosen so shared-category folders
# (multiple atlases living in one directory, e.g. YeoLab/WashU) get exercised too.
SPACE_CASES = [
    ("AL20", "LairdColin2mm"),
    ("TY17", "fsaverage6"),
    ("EG17", "fs_LR_32k"),
    ("EG5", "FSLMNI2mm"),
    ("XS268_8", "ShenColin1mm"),
]


def run_check(results, label, fn):
    """Run fn(), print PASS/FAIL, record the outcome in results. Never raises."""
    print(f"\n=== {label} ===")
    try:
        fn()
        print("  PASS")
        results.append((label, None))
    except Exception as exc:  # noqa: BLE001
        print(f"  FAIL: {exc!r}")
        traceback.print_exc()
        results.append((label, exc))


def check_atlas_fetch_per_space(cnc):
    for atlas_name, space in SPACE_CASES:
        params = cnc.compute_overlap_with_atlases.RefAtlasParams(atlas_name)
        assert params.config.space == space, (atlas_name, params.config.space, space)
        data = cnc.compute_overlap_with_atlases.read_atlas(params)
        assert data is not None and len(data) > 0, f"{atlas_name}: empty data"
        print(f"  {atlas_name} ({space}): OK")


def check_bad_atlas_name_raises(cnc):
    assert cnc.compute_overlap_with_atlases.if_atlas_exist("TY17") is True
    assert cnc.compute_overlap_with_atlases.if_atlas_exist("NOT_A_REAL_ATLAS") is False
    try:
        cnc.compute_overlap_with_atlases.compute_overlap_data_all(
            cnc.compute_overlap_with_atlases.RefAtlasParams("TY17"),
            ["NOT_A_REAL_ATLAS"],
        )
        raise AssertionError("expected ValueError for a nonexistent atlas")
    except ValueError:
        pass  # expected -- confirms the old sys.exit(1)/NameError bug is gone


def check_quick_example(cnc, out_dir):
    atlas_names_list = ["EG17", "MG360J12"]
    example = cnc.load_example
    ref_params = cnc.compute_overlap_with_atlases.DataParams(example.example_config, example.example_nii)
    result_dir = out_dir / "example_results"
    cnc.compute_overlap_with_atlases.network_correspondence(ref_params, atlas_names_list, str(result_dir))
    produced = list(result_dir.rglob("*"))
    assert produced, "network_correspondence produced no output files"
    print(f"  wrote {len(produced)} file(s) to {result_dir}")


def check_atlas_vs_atlas(cnc):
    overlap = cnc.compute_overlap_with_atlases.obtain_overlap_atlases("TY17", "EG17", "DefaultA", "Default")
    assert overlap is not None
    print(f"  TY17 vs EG17 overlap map: shape={overlap.shape}")


def check_load_and_draw_overlap_results(cnc, out_dir):
    overlap_val = cnc.compute_overlap_with_atlases.load_overlap_atlases("TY17", "EG17")
    assert overlap_val.size > 0
    figfile = out_dir / "overlap_mat_TY17_EG17.png"
    cnc.visualize_overlap_lib.draw_overlap_atlases("TY17", "EG17", 0, 1, str(figfile))
    assert figfile.exists(), f"figure was not written to {figfile}"
    print(f"  overlap matrix shape={overlap_val.shape}, figure -> {figfile}")


def check_spin_permutation(cnc):
    ref_params = cnc.compute_overlap_with_atlases.RefAtlasParams("TY17")
    p_val = cnc.compute_overlap_with_atlases.permute_overlap_data_all(ref_params, ["EG17"])
    assert "EG17" in p_val
    p = p_val["EG17"]
    assert ((p >= 0) & (p <= 1)).all(), "p-values out of [0, 1] range"
    print(f"  p-value matrix shape={p.shape}, range=[{p.min():.4f}, {p.max():.4f}]")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-spin", action="store_true",
                         help="skip the spin-rotation fetch/permutation test (saves ~1.4GB download)")
    args = parser.parse_args()

    t0 = time.time()
    import cbig_network_correspondence as cnc
    import_time = time.time() - t0
    print(f"import cbig_network_correspondence took {import_time:.2f}s")
    if import_time > 5:
        print("  NOTE: import took unexpectedly long -- it should do no network I/O by itself.")

    # sys.path ordering alone doesn't help if some earlier import already cached a
    # different cbig_network_correspondence in sys.modules (Python checks the
    # cache before consulting sys.path) -- so verify explicitly rather than assume.
    if not cnc.__file__.startswith(TOOLBOX_SRC):
        print(f"verify_release.py: imported cbig_network_correspondence from {cnc.__file__!r}, "
              f"not this local checkout ({TOOLBOX_SRC!r}) -- a different install is shadowing "
              "it. Start a fresh interpreter, or pip uninstall the released package from this "
              "environment if present.", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(tempfile.mkdtemp(prefix="cnc_verify_"))
    results = []

    run_check(results, "atlas fetch + read, one per space", lambda: check_atlas_fetch_per_space(cnc))
    run_check(results, "if_atlas_exist true/false + bad-name error path", lambda: check_bad_atlas_name_raises(cnc))
    run_check(results, "Quick Example from README.md", lambda: check_quick_example(cnc, out_dir))
    run_check(results, "atlas-vs-atlas comparison (cross-space)", lambda: check_atlas_vs_atlas(cnc))
    run_check(results, "load + draw precomputed overlap_results", lambda: check_load_and_draw_overlap_results(cnc, out_dir))

    if args.skip_spin:
        print("\n=== spin-rotation permutation test ===\n  SKIPPED (--skip-spin)")
    else:
        run_check(results, "spin-rotation permutation p-value (downloads ~1.4GB on first run)",
                   lambda: check_spin_permutation(cnc))

    shutil.rmtree(out_dir, ignore_errors=True)

    failures = [(label, exc) for label, exc in results if exc is not None]
    print("\n" + "=" * 50)
    if failures:
        print(f"{len(failures)} FAILURE(S):")
        for label, exc in failures:
            print(f"  - {label}: {exc!r}")
        sys.exit(1)
    print(f"All {len(results)} checks passed.")
    sys.exit(0)


if __name__ == "__main__":
    main()
