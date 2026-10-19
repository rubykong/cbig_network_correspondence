"""
Promote a staged atlas release into cbig_network_correspondence_data.

Run this only after inspecting the output of stage_atlas.sh and deciding it's
correct. This script:
  1. copies each staged atlas file into atlases/<space>/<category>/
  2. writes atlas_config/<Data_Name> from the config's [data_info] section, filling in
     Data_Threshold from [release_info] Thresholds if the user didn't set it directly
  3. appends Data_Name / Data_FullName to atlas_list / atlas_list_fullname (kept
     line-aligned -- both files are parallel arrays)
  4. copies Network_Names_File -> network_names/<Data_Name>, if given
  5. copies Network_Assignment_File -> network_assignment/<Data_Name>.mat, if given

It refuses to overwrite any existing file unless --force is given, and never runs
git itself -- add/commit/tag/push in cbig_network_correspondence_data is still a
manual step for you. overlap_results computation against existing atlases is
intentionally not part of this script (see README.md).

Usage: promote_atlas.py <config.ini> [staging_root] [--force]
"""
import argparse
import configparser
import shutil
import sys
from pathlib import Path

from parse_atlas_config import ReleaseConfig, _load_grab_data_info

DATA_REPO = Path("/home/rkong/storage/ruby/data/cbig_network_correspondence_data")
ALL_SPACES = ("fsaverage6", "fs_LR_32k", "FSLMNI2mm", "LairdColin2mm", "ShenColin1mm")


def copy_file(src, dst, force):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and not force:
        raise FileExistsError(f"{dst} already exists (use --force to overwrite)")
    shutil.copy2(src, dst)
    print(f"  wrote {dst}")


def promote_atlas_files(cfg, staging_root, force):
    copied_any = False
    for space in ALL_SPACES:
        staged_dir = staging_root / "atlases" / space / cfg.category
        if not staged_dir.is_dir():
            continue
        staged_files = [f for f in staged_dir.iterdir() if f.is_file()]
        matching = [f for f in staged_files if cfg.name in f.name]
        if not matching:
            continue
        dest_dir = DATA_REPO / "atlases" / space / cfg.category
        for staged_file in matching:
            copy_file(staged_file, dest_dir / staged_file.name, force)
            copied_any = True
    if not copied_any:
        raise FileNotFoundError(
            f"No staged files found for atlas {cfg.name!r} under {staging_root}/atlases/. "
            "Did stage_atlas.sh run (and succeed) for this config?"
        )


def _resolve_data_threshold(cfg, data_info):
    """Reconcile [release_info] Thresholds (already validated by ReleaseConfig to be
    uniform across files, see parse_atlas_config.py) with any Data_Threshold the user
    wrote directly into [data_info]. Returns the Data_Threshold string to write, or
    None if neither source gave one.

    Thresholds itself is never shipped into atlas_config -- it's a release-pipeline
    convenience for authoring a per-file-validated-uniform value; Data_Threshold is
    the one the toolbox actually reads at runtime (see grab_data_info.py)."""
    existing = data_info.get("Data_Threshold")
    if cfg.uniform_threshold is None:
        return existing
    if existing is None:
        return cfg.uniform_threshold
    grab_info = _load_grab_data_info()
    if grab_info.thresh_str_to_num(existing) != grab_info.thresh_str_to_num(cfg.uniform_threshold):
        raise ValueError(
            f"{cfg.name}: [release_info] Thresholds (-> {cfg.uniform_threshold}) and "
            f"[data_info] Data_Threshold ({existing}) disagree -- set only one of them"
        )
    return existing


def write_atlas_config(cfg, config_path, force):
    """Copy the [data_info] section verbatim -- not the grab_data_info-parsed and
    possibly-defaulted version -- so optional fields the user omitted stay omitted,
    except Data_Threshold is filled in from [release_info] Thresholds when the user
    didn't already set it themselves (see _resolve_data_threshold)."""
    parser = configparser.ConfigParser()
    parser.optionxform = str  # preserve key case (e.g. "Data_Category"); default lowercases
    parser.read(config_path)
    data_info = parser["data_info"]

    dest = DATA_REPO / "atlas_config" / cfg.name
    if dest.exists() and not force:
        raise FileExistsError(f"{dest} already exists (use --force to overwrite)")

    data_threshold = _resolve_data_threshold(cfg, data_info)
    lines = ["[data_info]"]
    wrote_threshold = False
    for key, value in data_info.items():
        lines.append(f"{key}: {value}")
        wrote_threshold = wrote_threshold or key == "Data_Threshold"
    if data_threshold is not None and not wrote_threshold:
        lines.append(f"Data_Threshold: {data_threshold}")
    dest.write_text("\n".join(lines) + "\n")
    print(f"  wrote {dest}")


def append_to_atlas_lists(cfg, force):
    list_path = DATA_REPO / "atlas_list"
    fullname_path = DATA_REPO / "atlas_list_fullname"

    names = list_path.read_text().splitlines()
    fullnames = fullname_path.read_text().splitlines()
    if len(names) != len(fullnames):
        raise ValueError(
            f"{list_path} and {fullname_path} are out of sync "
            f"({len(names)} vs {len(fullnames)} lines) -- fix that by hand before promoting"
        )

    if cfg.name in names:
        if not force:
            raise FileExistsError(
                f"{cfg.name!r} is already in atlas_list (use --force to re-append; "
                "this will NOT remove the earlier entry)"
            )
        print(f"  NOTE: {cfg.name!r} already in atlas_list, appending anyway (--force)")

    if not cfg.full_name:
        raise ValueError(
            f"{cfg.name}: Data_FullName is required in [release_info] to keep "
            "atlas_list_fullname aligned with atlas_list"
        )

    names.append(cfg.name)
    fullnames.append(cfg.full_name)
    list_path.write_text("\n".join(names) + "\n")
    fullname_path.write_text("\n".join(fullnames) + "\n")
    print(f"  appended {cfg.name!r} to atlas_list and {cfg.full_name!r} to atlas_list_fullname")


def write_network_names(cfg, force):
    if not cfg.network_names_file:
        return
    dest = DATA_REPO / "network_names" / cfg.name
    copy_file(Path(cfg.network_names_file), dest, force)


def write_network_assignment(cfg, force):
    if not cfg.network_assignment_file:
        return
    dest = DATA_REPO / "network_assignment" / f"{cfg.name}.mat"
    copy_file(Path(cfg.network_assignment_file), dest, force)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="path to the atlas release config.ini")
    parser.add_argument("staging_root", nargs="?", default=None,
                         help="defaults to Collaborators/NetworkName/staged_releases/<Data_Name>")
    parser.add_argument("--force", action="store_true",
                         help="overwrite existing files / re-append an already-listed atlas")
    args = parser.parse_args()

    try:
        cfg = ReleaseConfig(args.config)
        staging_root = Path(args.staging_root) if args.staging_root else Path(
            f"/home/rkong/storage/ruby/data/Collaborators/NetworkName/staged_releases/{cfg.name}"
        )
        if not staging_root.is_dir():
            raise FileNotFoundError(f"Staging root not found: {staging_root}. Run stage_atlas.sh first.")

        print(f"Promoting {cfg.name!r} from {staging_root} into {DATA_REPO}")

        print("Copying atlas files...")
        promote_atlas_files(cfg, staging_root, args.force)

        print("Writing atlas_config...")
        write_atlas_config(cfg, args.config, args.force)

        print("Updating atlas_list / atlas_list_fullname...")
        append_to_atlas_lists(cfg, args.force)

        print("Writing network_names / network_assignment (if provided)...")
        write_network_names(cfg, args.force)
        write_network_assignment(cfg, args.force)
    except (ValueError, FileNotFoundError, FileExistsError) as exc:
        print(f"promote_atlas.py: {exc}", file=sys.stderr)
        sys.exit(1)

    print()
    print(f"Done. {cfg.name!r} is now staged in your local {DATA_REPO} working copy.")
    print("This script did NOT run git. Review the changes, then from that directory:")
    print("  git status && git add -A && git commit -m 'Add atlas <name>' && git push")
    print("(and cut a new tag / bump DATA_TAG in the toolbox's datasets.py if you want")
    print(" the toolbox's lazy fetch to pick this up -- see scripts/README.md there.)")
    print()
    print("Not done by this script (see README.md): computing overlap_results against")
    print("existing atlases.")


if __name__ == "__main__":
    main()
