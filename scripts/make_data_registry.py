"""Generate Pooch registry files from a tagged commit of cbig_network_correspondence_data.

Usage:
    python scripts/make_data_registry.py <path-to-data-repo> <git-ref>

Reads the file list and content straight out of git (not the working tree),
so the resulting hashes are guaranteed to match whatever raw.githubusercontent.com
serves for that ref. Writes two registry files, each line `<relative-path> <sha256>`:
    src/cbig_network_correspondence/registry/core.txt            (everything)
    src/cbig_network_correspondence/registry/spin_rotations.txt  (spin_rotations/** only)
"""
import hashlib
import subprocess
import sys
from pathlib import Path

TOOLBOX_ROOT = Path(__file__).resolve().parent.parent
REGISTRY_DIR = TOOLBOX_ROOT / "src" / "cbig_network_correspondence" / "registry"


def git_file_list(repo, ref):
    out = subprocess.run(
        ["git", "-C", str(repo), "ls-tree", "-r", "--name-only", ref],
        check=True, capture_output=True, text=True,
    )
    return [line for line in out.stdout.splitlines() if line]


def git_file_hash(repo, ref, relpath):
    proc = subprocess.run(
        ["git", "-C", str(repo), "show", f"{ref}:{relpath}"],
        check=True, capture_output=True,
    )
    return hashlib.sha256(proc.stdout).hexdigest()


def main(repo, ref):
    repo = Path(repo)
    files = git_file_list(repo, ref)

    core_lines = []
    spin_lines = []
    for relpath in sorted(files):
        digest = git_file_hash(repo, ref, relpath)
        line = f"{relpath} {digest}\n"
        if relpath.startswith("spin_rotations/"):
            spin_lines.append(line)
        else:
            core_lines.append(line)

    REGISTRY_DIR.mkdir(parents=True, exist_ok=True)
    (REGISTRY_DIR / "core.txt").write_text("".join(core_lines))
    (REGISTRY_DIR / "spin_rotations.txt").write_text("".join(spin_lines))

    print(f"core.txt: {len(core_lines)} files")
    print(f"spin_rotations.txt: {len(spin_lines)} files")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
