# Releasing a new package version

There are two independent kinds of "release" for this project -- do them in this
order, since a package release ships whatever data tag is currently wired into
`datasets.py`:

1. **Data release** -- adding/updating an atlas in `cbig_network_correspondence_data`.
   See `scripts/release_atlas/README.md` for the full atlas-release workflow
   (`stage_atlas.sh` / `promote_atlas.py` / `compute_overlap.py`), then:
   ```bash
   python scripts/make_data_registry.py /path/to/cbig_network_correspondence_data vX.Y.Z
   ```
   and bump `DATA_TAG` in `src/cbig_network_correspondence/datasets.py` to match. Only
   needed when the data repo itself changed; skip straight to step 2 for a
   code-only release.

2. **Package release** -- cutting a new PyPI version of this toolbox (below).

Nothing in step 2 needs to be run locally -- `poetry` isn't even installed in most
dev environments here. The whole version-bump / changelog / tag / publish sequence
runs through GitHub Actions, triggered from the GitHub web UI.

## Before cutting a release

1. **Make sure `CHANGELOG.md`'s `[Unreleased]` section is accurate.** Every PR
   should already be adding its own entry there (Keep a Changelog format --
   `Added`/`Changed`/`Fixed`/etc.), but double check nothing is missing before you
   cut a version, since whatever's in `[Unreleased]` becomes the release notes.
2. **Run `scripts/verify_release.py`** in an environment with the full dependency
   stack (e.g. `NCT_env`/`Ruby_py3`):
   ```bash
   python scripts/verify_release.py              # full run, incl. spin-rotation test (~1.4GB)
   python scripts/verify_release.py --skip-spin   # quicker, skips that download
   ```
   See `scripts/README.md` for what each check covers. Fix anything that fails
   before continuing -- this is the last local gate before anything gets tagged
   or published.
3. If you changed anything touching packaging (`pyproject.toml`'s `include`,
   `datasets.py`'s registry loading, etc.), also do the wheel-packaging sanity
   check described at the bottom of `scripts/README.md` (`poetry build` +
   install the wheel into a throwaway venv + confirm `datasets.fetch_file()`
   resolves).

## Cutting the release (GitHub Actions)

1. On GitHub: **Actions -> "Draft a release" -> Run workflow**. Give it a
   `version` input -- either a bump keyword (`patch` / `minor` / `major` /
   `prepatch` / `preminor` / `premajor` / `prerelease`) or an explicit version
   like `1.2.3`; this is passed straight to `poetry version`.

   This workflow (`.github/workflows/draft_release.yml`) automatically:
   - bumps the version in `pyproject.toml` (`poetry version`)
   - moves `CHANGELOG.md`'s `[Unreleased]` section into a new dated
     `[X.Y.Z] - YYYY-MM-DD` section (`kacl-cli release --modify --auto-link`) --
     this is the "changelog thing" that tracks versions
   - commits `CHANGELOG.md` + `pyproject.toml` as `Release X.Y.Z`
   - creates and pushes git tag `X.Y.Z`
   - opens a **draft** GitHub release with that changelog section as the body

2. **Review the draft release** on GitHub (Releases page) -- check the version
   number and changelog body look right.
3. **Click "Publish release."** Publishing (not drafting) is what triggers
   `.github/workflows/release.yml`, which:
   - builds and publishes the package to PyPI (`poetry publish --build`, using
     the `PYPI_TOKEN` secret)
   - deploys the docs site (`mkdocs gh-deploy`)

That's it -- there's no local `poetry publish`/`twine upload` step; publishing
only happens by clicking "Publish release" on the draft GitHub Actions created.

## After publishing

Sanity-check the real PyPI artifact once it's up:
```bash
python -m venv /tmp/cnc_release_test
/tmp/cnc_release_test/bin/pip install cbig_network_correspondence
/tmp/cnc_release_test/bin/python -c "import cbig_network_correspondence as cnc; print(cnc.datasets.fetch_file('atlas_list'))"
```
If that prints a real local path (not an error), the published wheel picked up
the registry files correctly and can reach the tagged data.
