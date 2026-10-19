# Scripts

See `UPDATE_PACKAGE_README.md` for the full release checklist (data release vs.
package release, CHANGELOG/version bump, the GitHub Actions publish flow).

## `verify_release.py` — run this before cutting a PyPI release

Exercises the real package against the live, tagged `cbig_network_correspondence_data`
repo. This is the thing to run after making any change to the data-fetch layer
(`datasets.py`) or to `compute_overlap_with_atlases.py`/`visualize_overlap_lib.py`,
before publishing a new version.

### Prerequisites

This needs the *full* dependency set (`brainspace`, `nibabel`, `vtk`, etc.), so run it
inside the environment you actually use for this toolbox (e.g. `Ruby_py3` / `NCT_env`),
not a bare Python install. The dependencies need to already be installed there, but
`cbig_network_correspondence` itself does not need to be `pip install`-ed at all: the
script prepends this checked-out repo's `src/` to `sys.path` itself and verifies at
runtime that `cbig_network_correspondence.__file__` actually resolves there — so it
always tests your local changes, even if a released version is also installed
somewhere on the same environment's path.

You'll also need internet access — the first run downloads data from
`raw.githubusercontent.com/rubykong/cbig_network_correspondence_data` into a local
cache (`~/.cache/cbig_network_correspondence` by default, or
`$CBIG_NETWORK_CORRESPONDENCE_DATA_DIR` if you've set that). Subsequent runs reuse the
cache and are fast.

### Running it

```bash
python scripts/verify_release.py
```

This downloads the spin-rotation permutation data (~1.4GB total, across both surface
spaces) the first time it runs. If you just want a quick check and don't need to
re-verify the permutation-testing path, skip it:

```bash
python scripts/verify_release.py --skip-spin
```

### What it checks

1. **Import timing** — `import cbig_network_correspondence` should be near-instant
   (no network I/O at import time; data is only fetched when you actually use it).
2. **One atlas per space** (`AL20`/LairdColin2mm, `TY17`/fsaverage6, `EG17`/fs_LR_32k,
   `EG5`/FSLMNI2mm, `XS268_8`/ShenColin1mm) — fetches the atlas config, the atlas's
   data file(s), and reads them. Three of these five live in a folder shared with
   other atlases (`YeoLab`, `WashU`), which exercises the name-filtering in
   `datasets.fetch_dir`.
3. **Error path for a nonexistent atlas** — confirms a bad atlas name raises
   `ValueError` (not the old `NameError` bug from the missing `sys` import).
4. **The Quick Example from `README.md`** — the exact snippet documented for new
   users, run for real, checking it actually produces output files.
5. **An atlas-vs-atlas comparison** across two different atlases/categories.
6. **Loading + drawing a precomputed `overlap_results` entry** — exercises
   `visualize_overlap_lib.draw_overlap_atlases` and confirms a figure file gets written.
7. **Spin-rotation permutation testing** (unless `--skip-spin`) — fetches and
   reassembles the part-split `spin_rotations` data for a space, runs a real
   permutation test, and checks the resulting p-values fall in `[0, 1]`.

### Reading the output

The script prints `PASS`/`FAIL` per check as it goes, a full traceback for anything
that fails, and a summary at the end. It exits with code `0` if everything passed and
`1` otherwise — safe to use as a gate (e.g. `python scripts/verify_release.py && echo
"ok to release"`).

### One more check this script can't do: the actual release artifact

`pip install -e .` reads files straight from this source tree, so it will not catch a
packaging mistake — specifically, whether the `registry/core.txt` and
`registry/spin_rotations.txt` files (which `datasets.py` needs to find anything at all)
actually get bundled into a real build. Before publishing, also do this once:

```bash
poetry build
python -m zipfile -l dist/*.whl | grep registry   # should list both registry/*.txt files
python -m venv /tmp/cnc_pypi_test
/tmp/cnc_pypi_test/bin/pip install dist/*.whl
/tmp/cnc_pypi_test/bin/python -c "import cbig_network_correspondence as cnc; print(cnc.datasets.fetch_file('atlas_list'))"
```

If that last line prints a real local path instead of erroring, the wheel is correctly
packaged for release.

## `make_data_registry.py`

Regenerates `src/cbig_network_correspondence/registry/{core,spin_rotations}.txt` from
a git ref of the `cbig_network_correspondence_data` repo. Only needs to be re-run when
the data repo gets a new tagged release (and `datasets.DATA_TAG` is bumped to match):

```bash
python scripts/make_data_registry.py /path/to/cbig_network_correspondence_data v1.1.0
```

It reads file contents straight from git (`git show <ref>:<path>`), not the working
tree, so the hashes it writes are guaranteed to match whatever
`raw.githubusercontent.com` serves for that ref regardless of local working-tree state.
