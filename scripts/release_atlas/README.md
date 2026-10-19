# Atlas release pipeline

Replaces the old workflow of copy-pasting one of the per-atlas scripts under
`Collaborators/NetworkName/scripts/` and hand-editing filenames inside it. Now:
one config file describes the atlas, `stage_atlas.sh` projects it into every
supported space under a staging directory you can inspect, `promote_atlas.py`
copies the result into `cbig_network_correspondence_data` once you're satisfied,
and `compute_overlap.py` fills in `overlap_results` against the atlases already
there.

```
stage_atlas.sh <config.ini> [staging_root]      # project, write to staging only
  ... inspect the staged files ...
promote_atlas.py <config.ini> [staging_root]    # copy staging -> the real data repo
  ... check the promoted atlas looks right (see "After promoting" below) ...
compute_overlap.py <config.ini>                 # fill in overlap_results vs existing atlases
  ... git add/commit/tag/push yourself when ready ...
```

None of the three scripts run git. Each only touches your local
`cbig_network_correspondence_data` working copy.

## Supported native spaces

| Native space      | Atlas type supported | What happens |
|---|---|---|
| `FSLMNI2mm`        | Hard / Soft / Metric | Generic MNI fan-out directly (`proj_mni_fanout.m`) |
| `fs_LR_32k`        | Hard only            | Surface projection to the other surface space + MNI (`proj_native_to_mni.m`), then volume fan-out |
| `fsaverage6`       | Hard only            | Same as above, with an extra fsaverage6→fsaverage7 upsample step |
| `LairdColin2mm`    | Hard / Soft / Metric | Forward ANTs warp to MNI (`proj_volume_native_to_mni.sh`), then the generic MNI fan-out |
| `ShenColin1mm`     | Hard / Soft / Metric | Same as above, with Shen's own warp |
| anything else      | --                    | `stage_atlas.sh` refuses with a clear error before running anything |

fs_LR_32k/fsaverage6-native is Hard-only because the two scripts it generalizes
(`proj_Glasser_to_fs_MNI.m`, `proj_Yeo_to_fslr_MNI.m`) only ever handled a single
`lh_labels`/`rh_labels` map. If you have a Soft/Metric atlas that's natively a
surface map, convert it to an `FSLMNI2mm` `.nii.gz` yourself first (exactly as was
already being done for atlases like HCPICA) and release it as FSLMNI2mm-native.

**Not automated, still your job**: converting a raw source file (cifti, a
collaborator's custom format, whatever) into the canonical native-space input this
pipeline expects (a `.nii.gz` for volumetric native spaces, or a `.mat` with
`lh_labels`/`rh_labels` for surface native spaces). That conversion varies too much
per data source to generalize -- see the existing one-off `write_*_as_mat.m`
scripts for examples of what that conversion has looked like before.

## Config file format

One INI file per atlas. See `atlas_release_config.example.ini` for an annotated
template. Two sections:

### `[data_info]`

Identical format to the files already in `cbig_network_correspondence_data/atlas_config/`
-- in fact `promote_atlas.py` copies this section verbatim into
`atlas_config/<Data_Name>`. Use the same casing (`Data_Category`, not
`data_category`) and **omit a key entirely rather than leaving it blank** -- e.g.
don't write `Data_Threshold:` with nothing after the colon for a Hard atlas; that
crashes `grab_data_info`'s parser (it tries to split an empty string into two
threshold bounds). Look at an existing `atlas_config` file matching your `Data_Type`
for a template.

### `[release_info]`

| Field | Required | Meaning |
|---|---|---|
| `Data_FullName` | yes | Human-readable name, appended to `atlas_list_fullname` alongside `Data_Name` in `atlas_list` (these two files are parallel arrays and must stay line-aligned) |
| `Source_Files` | yes | Path or glob to the native-space file(s). For `fs_LR_32k`/`fsaverage6` this must match exactly one file. Each matched filename must contain `Data_Name` as a substring -- that's how the toolbox finds an atlas's file(s) at runtime (`sort_files`'s substring match). |
| `Thresholds` | no | Semicolon-separated threshold strings (e.g. `[1.5,Inf]; [1.5,Inf]`), one per file matched by `Source_Files`, **in natural-sorted (natsort) order**. Only valid for `FSLMNI2mm`/`LairdColin2mm`/`ShenColin1mm`-native. **Not applied during staging** -- staged/promoted atlas data stays exactly as continuous as the source files. All entries must parse to the same `(low, high)` pair, since the toolbox's `Data_Threshold` (what this actually feeds) is a single value applied uniformly to every file of a Soft/Metric atlas at runtime, not a per-file one -- `parse_atlas_config.py` errors if your entries disagree. `promote_atlas.py` writes that single value into `Data_Threshold` in the promoted `atlas_config` (unless you already set `Data_Threshold` yourself in `[data_info]`, in which case the two must agree). |
| `Network_Names_File` | yes | One network name per line (same order as the parcellation labels, or the natural-sorted `Source_Files` for a multi-map Soft/Metric atlas -- e.g. an ICA decomposition with one file per component), copied to `network_names/<Data_Name>`. Every atlas currently in `cbig_network_correspondence_data` has one -- it's what `construct_network_name()` reads at runtime. For Soft/Metric atlases, the line count must match the number of files matched by `Source_Files`. |
| `Network_Assignment_File` | no | A `.mat` with a `mapping` variable, copied to `network_assignment/<Data_Name>.mat`. Required if `Data_NetworkAssignment` is true. |

`parse_atlas_config.py` validates all of this up front -- unsupported native space,
wrong file count for a surface-native atlas, `Thresholds` count not matching
`Source_Files` count, missing `Network_Assignment_File` when
`Data_NetworkAssignment` is true, etc. -- and fails with a specific message before
`stage_atlas.sh` runs anything expensive.

**Ordering**: `Thresholds` and `Network_Names_File` are matched to `Source_Files` by
**natural sort** (`parse_atlas_config.py` bundles its own stdlib-only natural-sort
function, so this has no extra dependency to install) -- matching the order
`sort_files()`'s `natsorted()` (from the `natsort` package) uses in
`compute_overlap_with_atlases.py`.
That means `zstat1, zstat2, ..., zstat9, zstat10, zstat11, zstat12`, not the
lexicographic `zstat1, zstat10, zstat11, zstat12, zstat2, ...` that plain `sorted()`
or a naive `dir()` glob would produce. This matters because that's exactly how most
ICA decompositions get named (one unpadded-numbered file per component), and the
toolbox will read these files back in natural-sorted order at runtime too -- so
getting this wrong wouldn't error, it would silently pair the wrong network name
with the wrong map. `stage_atlas.sh` resolves the natural-sorted file list once (in
Python) and passes that exact list into MATLAB explicitly, rather than letting
MATLAB's `dir()` re-resolve and potentially re-sort it independently. If you're ever
unsure what order your files will be matched in, check with
`python3 parse_atlas_config.py your_config.ini` before staging -- it prints the
resolved file list in the order `Thresholds`/`Network_Names_File` will be matched
against.

## Worked examples

### FSLMNI2mm-native (any type)

```ini
[data_info]
Data_Category: Gordon
Data_Name: WashU_12Networks
Data_Space: FSLMNI2mm
Data_Type: Hard

[release_info]
Data_FullName: Gordon/WashU 12-network parcellation
Source_Files: /home/rkong/storage/ruby/data/Collaborators/NetworkName/Gordon/MNI/WashU_12Networks.nii.gz
Network_Names_File: /home/rkong/storage/ruby/data/Collaborators/NetworkName/Gordon/WashU_12Networks_network_names.txt
```

```bash
bash stage_atlas.sh washu_12networks.ini
# inspect staging_root/atlases/{fsaverage6,fs_LR_32k,FSLMNI2mm,LairdColin2mm,ShenColin1mm}/Gordon/
python3 promote_atlas.py washu_12networks.ini
```

### FSLMNI2mm-native, one file per network (the ICA pattern)

The common case for an ICA decomposition, dual regression maps, or any other
atlas that's naturally "one map per network" rather than one multi-label file:
`Source_Files` is a glob matching every component map, and `Thresholds`/
`Network_Names_File` each have one entry per matched file. This is Soft/Metric,
since each file is its own continuous spatial map, not a single hard-labeled
parcellation. The staged/promoted data stays exactly as continuous as these
source files -- `Thresholds` here (all `[1.5,Inf]`, i.e. uniform across every
component) only sets `Data_Threshold` in the promoted `atlas_config`, applied by
the toolbox at *read* time (e.g. when computing `overlap_results`).

```ini
[data_info]
Data_Category: HCPICA
Data_Name: HCPICA_d25
Data_Space: FSLMNI2mm
Data_Type: Soft

[release_info]
Data_FullName: HCP ICA, 25-component decomposition
Source_Files: /path/to/native/HCPICA_d25_thresh_zstat*.nii.gz
Thresholds: [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; [1.5,Inf]
Network_Names_File: /path/to/HCPICA_d25_network_names.txt   # 25 lines, one per component, natural-sorted order
```

Filenames like `HCPICA_d25_thresh_zstat1.nii.gz` ... `HCPICA_d25_thresh_zstat25.nii.gz`
are naturally sorted (component 1, 2, 3, ..., 25), *not* alphabetically (which would
put component 10 right after component 1) -- see "Ordering" above. Run
`python3 parse_atlas_config.py hcpica.ini` first if you want to double-check the
exact file order your `Thresholds`/`Network_Names_File` entries need to line up
with before staging.

The same pattern applies to a multi-file atlas natively in `LairdColin2mm` or
`ShenColin1mm` -- see the Laird example below.

### fs_LR_32k-native, Hard only (the Glasser pattern)

```ini
[data_info]
Data_Category: Glasser
Data_Name: MG360J12
Data_Space: fs_LR_32k
Data_Type: Hard

[release_info]
Data_FullName: Glasser MMP 360-region, Ji 2019 12-network solution
Source_Files: /path/to/native_fs_LR_32k/MG360J12.mat
Network_Names_File: /path/to/MG360J12_network_names.txt
```

Same two commands as above. `stage_atlas.sh` writes the fsaverage6 `.mat` and an
MNI `.nii.gz` via `proj_native_to_mni.m`, then runs `proj_mni_fanout_volumes.sh`
directly for LairdColin2mm/ShenColin1mm -- it deliberately does *not* re-run the
full MATLAB fan-out for this case, since that would overwrite the
surface-to-surface fsaverage6 projection (more accurate for a genuinely
surface-native atlas) with a lower-fidelity MNI round-trip.

### fsaverage6-native, Hard only (the Yeo pattern)

Same as above with `Data_Space: fsaverage6`. Internally this takes one extra step
(upsampling fsaverage6 to full fsaverage7 resolution via nearest-neighbor
interpolation) before the MNI projection -- same approximation the original
`proj_Yeo_to_fslr_MNI.m` already made, not something new.

### LairdColin2mm-native or ShenColin1mm-native (new)

```ini
[data_info]
Data_Category: Laird
Data_Name: Laird_ICA
Data_Space: LairdColin2mm
Data_Type: Soft

[release_info]
Data_FullName: Laird 2011 ICA, 20 components
Source_Files: /path/to/native_LairdColin2mm/Laird_ICA_zstat*.nii.gz
Thresholds: [1.5,Inf]; [1.5,Inf]; [1.5,Inf]; ...   # one per matched file, sorted order
Network_Names_File: /path/to/Laird_ICA_network_names.txt   # 20 lines, one per zstat file, sorted order
```

`stage_atlas.sh` forward-warps each matched file to FSLMNI2mm
(`proj_volume_native_to_mni.sh`), then runs the generic MNI fan-out. The generic
fan-out also reconstructs a LairdColin2mm version via a volumetric round trip
(native → MNI → native) as a side effect of warping to *every* other space --
`stage_atlas.sh` automatically overwrites that with your real original
LairdColin2mm file(s) afterward, since a round trip is strictly lower fidelity than
the atlas you actually have. You'll see a
`Restored authentic native-space file(s)...` line confirming this.

## Recommended first real test

Before trusting this on a brand-new atlas, run it on one you've **already released
by hand** so you have known-good output to diff against:

```bash
bash stage_atlas.sh glasser_mg360j12.ini /tmp/stage_test
diff -rq /tmp/stage_test/atlases cbig_network_correspondence_data/atlases   # expect no diffs
# or, for .mat files where you want a numeric comparison rather than byte-equality:
#   compare lh_labels/rh_labels in MATLAB between the two
```

The config-parsing/dispatch logic (`parse_atlas_config.py`, `stage_atlas.sh`'s
branching, `promote_atlas.py`, `compute_overlap.py`) has since been tested for
real, including `compute_overlap.py`'s Dice computation against real atlas pairs
with known ground truth (see "After promoting" below). The actual MATLAB/ANTs
registration math in `proj_mni_fanout.m`, `proj_native_to_mni.m`, and the two
`.sh` warp scripts still has not been numerically verified end-to-end by me --
everything there was written by closely mirroring scripts that have already run
successfully in production (see the comments at the top of each file for which
original script it generalizes, or which existing forward-warp call it mirrors).
Please run the hand-release comparison above before trusting a brand-new atlas's
projected output.

## After promoting: checking the atlas, then computing overlap_results

Before running `compute_overlap.py`, check the promoted atlas is actually correct
-- a bad `atlas_config` won't necessarily error, it'll just silently produce wrong
results. In particular:

- **`Data_Type` must match the file layout.** `Hard` means a single multi-label
  file; `read_atlas()` only ever reads the *first* file `sort_files()` matches,
  even if `Source_Files` matched several. A multi-file atlas -- one file per
  network, e.g. an ICA decomposition or a set of binary RSN masks -- needs
  `Soft` (or `Metric`), not `Hard`, or everything past the first network is
  silently dropped. Sanity-check with:
  ```bash
  python3 -c "
  from cbig_network_correspondence import compute_overlap_with_atlases as coa
  print(len(coa.sort_files('cbig_network_correspondence_data/atlases/<space>/<category>', '<Data_Name>')), 'files matched')
  "
  ```
  and confirm that count matches your `Network_Names_File` line count (already
  enforced by `parse_atlas_config.py` for Soft/Metric, but worth eyeballing
  directly against the promoted files too).
- Spot-check one known comparison (e.g. `RefAtlasParams('<Data_Name>')` +
  `read_atlas()`, or an `obtain_overlap_atlases` call) the way
  `scripts/verify_release.py` does, before computing overlap against everything.

Once you're satisfied:

```bash
python3 compute_overlap.py my_atlas.ini                       # vs every atlas in atlas_list
python3 compute_overlap.py my_atlas.ini --atlases TY17,EG17    # vs specific atlases only
python3 compute_overlap.py my_atlas.ini --force                # recompute existing pairs too
```

This needs the full toolbox dependency stack (unlike `parse_atlas_config.py`/
`promote_atlas.py`), since it reuses the toolbox's own `compute_overlap()`/
`read_atlas()` for the actual Dice math -- run it in the same environment you use
for the toolbox itself. For each existing atlas it writes **both** directions,
`overlap_results/<new>/<existing>.mat` and `overlap_results/<existing>/<new>.mat`,
matching the convention every existing pair already follows -- just the Dice
matrix under the `overlap` key (no p-values; that's all `load_overlap_atlases`/
`draw_overlap_atlases` ever read). It skips a pair whose file already exists
unless `--force`, and operates on your local `cbig_network_correspondence_data`
checkout directly (not via `datasets.py`'s tagged-release fetch), since a
just-promoted atlas isn't tagged yet.

I validated this against real atlas pairs with pre-existing `overlap_results` as
ground truth (correlation ≥0.999 across every pair tested, exact min/max match on
several). Small numeric drift (up to ~0.09 in the volumetric-space pairs) showed
up against the stored files -- consistent with this repo's own history of
periodic "regenerate overlap mat" commits (e.g. after nifti header fixes), not a
bug in this script: it recomputes from current source data using the toolbox's
current code, same as those regeneration commits do.

### Checking the result visually

`Visualize_Overlap_Results.ipynb` (same directory) draws every
`overlap_results/<NEW_ATLAS>/*.mat` file on disk as a heatmap -- set `NEW_ATLAS`
to the atlas you just promoted and run it after `compute_overlap.py`, as a quick
visual sanity check (sane values, expected structure, no empty/garbage matrices)
before committing. Like `compute_overlap.py`, it's local-only (no network) and
always imports this checked-out repo's code, never a `pip install`-ed version.

## Explicitly out of scope

- **git operations.** None of the three scripts run `git add`/`commit`/`tag`/
  `push`. See "Git operations: getting lazy-fetch to serve the new atlas" below.
- **Raw-format-to-canonical-native conversion** (see above).

## Git operations: getting lazy-fetch to serve the new atlas

`datasets.py` in the toolbox repo fetches everything from
`raw.githubusercontent.com/rubykong/cbig_network_correspondence_data/<DATA_TAG>/`,
verifying each file against a hash recorded in
`src/cbig_network_correspondence/registry/{core,spin_rotations}.txt`. Promoting
an atlas with `promote_atlas.py`/`compute_overlap.py` only changes your local
`cbig_network_correspondence_data` working copy -- raw.githubusercontent won't
serve any of it, and the registry won't have hashes for it, until you do the
following, in order.

### 1. Push the data repo and cut a new tag

In `cbig_network_correspondence_data`:

```bash
git add atlas_config/<Data_Name> atlas_list atlas_list_fullname \
        network_names/<Data_Name> atlases/*/<Data_Category>/<Data_Name>* \
        overlap_results/<Data_Name>
# plus network_assignment/<Data_Name>.mat if you set Network_Assignment_File
git commit -m "Add <Data_Name>"
git push origin <your-branch>
# open a PR into master and merge it, same as existing atlases -- then on master:
git pull origin master
git tag v1.1.0 -m "Add <Data_Name>"       # bump the version, don't reuse v1.0.0
git push origin v1.1.0
```

The tag has to exist on GitHub (`git push origin <tag>`) before the next step --
`raw.githubusercontent.com` serves files by tag/ref, and `make_data_registry.py`
reads the tagged tree with `git show <ref>:<path>`.

### 2. Regenerate the registry and bump `DATA_TAG` in the toolbox repo

In `cbig_network_correspondence` (the toolbox repo):

```bash
python scripts/make_data_registry.py /path/to/cbig_network_correspondence_data v1.1.0
```

This overwrites `src/cbig_network_correspondence/registry/core.txt` and
`registry/spin_rotations.txt` with hashes for every file at that tag. Then bump
the version string in `src/cbig_network_correspondence/datasets.py`:

```python
DATA_TAG = "v1.1.0"
```

### 3. Verify, then commit/push the toolbox repo

```bash
python scripts/verify_release.py   # confirms the new tag's files actually fetch and read correctly
git add src/cbig_network_correspondence/registry/core.txt \
        src/cbig_network_correspondence/registry/spin_rotations.txt \
        src/cbig_network_correspondence/datasets.py
git commit -m "Bump DATA_TAG to v1.1.0 for <Data_Name>"
git push origin <your-branch>   # PR + merge as usual
```

At this point anyone running this toolbox repo's checked-out source (or a fresh
`pip install` once a new PyPI release is cut per `scripts/README.md`) will
lazy-fetch the new atlas on first use. Until `DATA_TAG` is bumped, `datasets.py`
keeps pointing at the old tag and the new atlas stays invisible to lazy-fetch
even though it's sitting in `cbig_network_correspondence_data` on GitHub.
