"""
Shared config reader for the atlas release pipeline (stage_atlas.sh / promote_atlas.py).

A release config is one INI file with two sections:

    [data_info]       -- identical format to the atlas_config/<name> files already
                          shipped in cbig_network_correspondence_data; this section is
                          copied verbatim into atlas_config/<Data_Name> at promote time.
    [release_info]     -- release-pipeline-only fields (not shipped as-is):
                             Data_FullName            -- human-readable full name,
                                                          appended to atlas_list_fullname
                             Source_Files             -- path or glob to the atlas's
                                                          native-space file(s)
                             Thresholds               -- optional, semicolon-separated
                                                          list of threshold strings (e.g.
                                                          '[1.5,Inf]; [1.5,Inf]'), one per
                                                          file matched by Source_Files, in
                                                          naturally sorted (natsort) order, matching sort_files() in
                                                          compute_overlap_with_atlases.py. Only used
                                                          for FSLMNI2mm/LairdColin2mm/
                                                          ShenColin1mm-native atlases.
                                                          Not applied at release time --
                                                          staged/promoted atlas data stays
                                                          continuous. All entries must parse
                                                          to the same (low, high) pair
                                                          (the toolbox's Data_Threshold is a
                                                          single global value applied to
                                                          every file of a Soft/Metric atlas,
                                                          see grab_data_info.py), and that
                                                          single value is written into
                                                          Data_Threshold in the promoted
                                                          atlas_config by promote_atlas.py,
                                                          so the toolbox applies it at
                                                          runtime (e.g. when computing
                                                          overlap_results).
                             Network_Names_File       -- required, one network name per
                                                          line (in the same order as the
                                                          parcellation labels / soft map
                                                          files), copied to network_names/
                                                          <name>. Every atlas currently in
                                                          cbig_network_correspondence_data
                                                          has one of these -- it's what
                                                          construct_network_name() reads.
                             Network_Assignment_File  -- optional .mat, copied to
                                                          network_assignment/<name>.mat
                                                          (only meaningful when
                                                          Data_NetworkAssignment is true)

See README.md in this directory for the full schema and a worked example.
"""
import configparser
import contextlib
import glob
import importlib.util
import io
import re
from pathlib import Path

_DIGIT_RUN = re.compile(r"(\d+)")


def _natsort_key(value):
    """Split into alternating text/number chunks so number chunks compare
    numerically (zstat2 < zstat10), not lexicographically (zstat10 < zstat2).
    Deliberately bundled rather than depending on the `natsort` PyPI package, so
    this script has no dependency beyond the standard library -- it only needs to
    match sort_files()'s natsorted() in compute_overlap_with_atlases.py closely
    enough for ordinary filenames, which this standard split-and-compare recipe does."""
    return [int(chunk) if chunk.isdigit() else chunk.lower() for chunk in _DIGIT_RUN.split(str(value))]


def natsorted(items):
    return sorted(items, key=_natsort_key)


SUPPORTED_SPACES = ("FSLMNI2mm", "fs_LR_32k", "fsaverage6", "LairdColin2mm", "ShenColin1mm")

# native spaces where the Glasser/Yeo surface-projection trick (proj_native_to_mni.m)
# applies -- it assumes a single lh/rh-label .mat, i.e. a Hard parcellation.
SURFACE_NATIVE_SPACES = ("fs_LR_32k", "fsaverage6")

# native spaces that need a forward ANTs warp to FSLMNI2mm before the MNI fan-out.
VOLUME_NATIVE_SPACES_NEEDING_WARP = ("LairdColin2mm", "ShenColin1mm")

TOOLBOX_SRC = Path(__file__).resolve().parents[2] / "src" / "cbig_network_correspondence"


def _load_grab_data_info():
    """
    Import grab_data_info.py directly by file path, bypassing
    cbig_network_correspondence/__init__.py (which would otherwise require
    brainspace/nibabel/vtk just to parse a config file).
    """
    spec = importlib.util.spec_from_file_location("grab_data_info", TOOLBOX_SRC / "grab_data_info.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReleaseConfig:
    def __init__(self, config_path):
        config_path = Path(config_path)
        if not config_path.is_file():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        self._grab_info = _load_grab_data_info()
        # read_config() prints informational notices (e.g. default-threshold notes) to
        # stdout; suppress them so stdout stays safe to `eval` (see _emit_shell below).
        with contextlib.redirect_stdout(io.StringIO()):
            self.data_info = self._grab_info.read_config(str(config_path))

        parser = configparser.ConfigParser()
        parser.read(config_path)
        if not parser.has_section("release_info"):
            raise ValueError(f"{config_path}: missing [release_info] section")
        release = parser["release_info"]
        self.full_name = release.get("Data_FullName", "").strip()
        self.source_files_pattern = release.get("Source_Files", "").strip()
        thresholds_raw = release.get("Thresholds", "").strip()
        self.thresholds = [t.strip() for t in thresholds_raw.split(";") if t.strip()] or None
        # set by _validate() once all Thresholds entries are confirmed to parse to the
        # same (low, high) pair -- the raw string of that one shared value, suitable
        # for writing straight into Data_Threshold in the promoted atlas_config.
        self.uniform_threshold = None
        self.network_names_file = release.get("Network_Names_File", "").strip() or None
        self.network_assignment_file = release.get("Network_Assignment_File", "").strip() or None

        self._validate()

    def _validate(self):
        name = self.data_info.name
        space = self.data_info.space
        dtype = self.data_info.type

        if not name:
            raise ValueError("Data_Name is required")
        if space not in SUPPORTED_SPACES:
            raise ValueError(
                f"Unsupported native space {space!r} for atlas {name!r}. "
                f"Supported native spaces: {', '.join(SUPPORTED_SPACES)}. "
                "If your atlas is natively in a different space, project it into one "
                "of these yourself first (see README.md)."
            )
        if space in SURFACE_NATIVE_SPACES and dtype != "Hard":
            raise ValueError(
                f"Native space {space!r} is only supported for Hard parcellations "
                f"(atlas {name!r} has Data_Type={dtype!r}). Soft/Metric atlases must "
                "be converted to an FSLMNI2mm .nii.gz yourself first -- see README.md."
            )
        if not self.source_files_pattern:
            raise ValueError(f"{name}: Source_Files is required in [release_info]")
        if space in SURFACE_NATIVE_SPACES:
            n_files = len(self.resolve_source_files())
            if n_files != 1:
                raise ValueError(
                    f"{name}: native space {space!r} expects exactly one Source_Files "
                    f"match (a single lh/rh-label .mat), got {n_files}"
                )
        if not self.network_names_file:
            raise ValueError(
                f"{name}: Network_Names_File is required in [release_info] -- every atlas "
                "currently in cbig_network_correspondence_data has one, since it's what "
                "construct_network_name() reads at runtime"
            )
        if not Path(self.network_names_file).is_file():
            raise FileNotFoundError(f"{name}: Network_Names_File not found: {self.network_names_file}")
        if self.data_info.netassign and not self.network_assignment_file:
            raise ValueError(
                f"{name}: Data_NetworkAssignment is true but no Network_Assignment_File given"
            )
        if self.network_assignment_file and not Path(self.network_assignment_file).is_file():
            raise FileNotFoundError(
                f"{name}: Network_Assignment_File not found: {self.network_assignment_file}"
            )
        if self.thresholds is not None:
            if space in SURFACE_NATIVE_SPACES:
                raise ValueError(
                    f"{name}: Thresholds is not applicable for surface-native space {space!r} "
                    "(only Hard parcellations are supported there, which aren't thresholded)"
                )
            n_files = len(self.resolve_source_files())
            if len(self.thresholds) != n_files:
                raise ValueError(
                    f"{name}: Thresholds has {len(self.thresholds)} entries but Source_Files "
                    f"matched {n_files} file(s) -- these must match 1:1, in natural-sorted order"
                )
            # Data_Threshold (the toolbox's runtime threshold, see grab_data_info.py's
            # thresh_str_to_num) is a single (low, high) pair applied uniformly to every
            # file of a Soft/Metric atlas -- it has no per-file variation. So a per-file
            # Thresholds list here only makes sense if every entry agrees; promote_atlas.py
            # writes that single value into Data_Threshold.
            parsed = [self._grab_info.thresh_str_to_num(t) for t in self.thresholds]
            if len(set(parsed)) > 1:
                raise ValueError(
                    f"{name}: Thresholds entries differ ({', '.join(self.thresholds)}) -- "
                    "Data_Threshold in [data_info] is a single value the toolbox applies to "
                    "every file of a Soft/Metric atlas at runtime (see grab_data_info.py's "
                    "thresh_str_to_num), so per-file thresholds can't be represented. Use the "
                    "same threshold for every file, or omit Thresholds here and set "
                    "Data_Threshold directly in [data_info] if you need something more exotic."
                )
            self.uniform_threshold = self.thresholds[0]
        if dtype in ("Soft", "Metric") and space not in SURFACE_NATIVE_SPACES:
            # one network name per native file for a multi-map soft/metric atlas,
            # matching the convention every existing Soft/Metric atlas already follows
            n_files = len(self.resolve_source_files())
            n_names = len(Path(self.network_names_file).read_text().splitlines())
            if n_names != n_files:
                raise ValueError(
                    f"{name}: Network_Names_File has {n_names} line(s) but Source_Files "
                    f"matched {n_files} file(s) -- for a Soft/Metric atlas these should "
                    "match 1:1"
                )

    def resolve_source_files(self):
        # natural sort, not lexicographic: must match sort_files()'s natsorted() in
        # compute_overlap_with_atlases.py exactly, since that's the order the toolbox
        # will later associate with Network_Names_File's lines at runtime. Lexicographic
        # sort would put "zstat10" before "zstat2" for an unpadded multi-file ICA-style
        # atlas, silently misaligning Thresholds/Network_Names_File against the wrong files.
        matches = natsorted(glob.glob(self.source_files_pattern))
        if not matches:
            raise FileNotFoundError(
                f"Source_Files pattern matched no files: {self.source_files_pattern!r}"
            )
        return matches

    @property
    def name(self):
        return self.data_info.name

    @property
    def category(self):
        return self.data_info.category

    @property
    def space(self):
        return self.data_info.space

    @property
    def dtype(self):
        return self.data_info.type


def _emit_shell(cfg):
    """Print `KEY='value'` lines, shell-quoted, for `eval` by stage_atlas.sh."""
    import shlex

    def line(key, value):
        print(f"{key}={shlex.quote('' if value is None else str(value))}")

    line("CNC_NAME", cfg.name)
    line("CNC_CATEGORY", cfg.category)
    line("CNC_SPACE", cfg.space)
    line("CNC_TYPE_LOWER", cfg.dtype.lower())
    line("CNC_SOURCE_FILES_PATTERN", cfg.source_files_pattern)
    # Natural-sorted, explicit file list -- bash array literal. stage_atlas.sh passes
    # this straight through to MATLAB instead of a glob, so the order is guaranteed
    # identical to the order validated here (and to what sort_files()'s natsorted()
    # will use when the toolbox reads this atlas back). Thresholds is deliberately
    # not emitted here -- stage_atlas.sh no longer applies it (see Thresholds'
    # docstring above); promote_atlas.py reads it from the ReleaseConfig object
    # directly instead.
    source_files = cfg.resolve_source_files()
    quoted = " ".join(shlex.quote(f) for f in source_files)
    print(f"CNC_SOURCE_FILES=({quoted})")


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 3 and sys.argv[1] == "--emit-shell":
        _emit_shell(ReleaseConfig(sys.argv[2]))
    elif len(sys.argv) == 2:
        cfg = ReleaseConfig(sys.argv[1])
        print(f"name={cfg.name} category={cfg.category} space={cfg.space} type={cfg.dtype}")
        print(f"source files: {cfg.resolve_source_files()}")
        print(f"thresholds: {cfg.thresholds} (uniform value for Data_Threshold: {cfg.uniform_threshold})")
    else:
        print(f"Usage: {sys.argv[0]} [--emit-shell] <config.ini>", file=sys.stderr)
        sys.exit(2)
