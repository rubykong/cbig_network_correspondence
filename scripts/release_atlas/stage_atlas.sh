#!/usr/bin/env bash
# Stage an atlas for release: projects it from its native space into every other
# supported space, writing everything under a staging directory that mirrors
# cbig_network_correspondence_data's layout -- nothing touches the real data repo
# here. Inspect the staged output, then run promote_atlas.py when you're satisfied.
#
# Usage: stage_atlas.sh <config.ini> [staging_root]
#   staging_root defaults to
#   /home/rkong/storage/ruby/data/Collaborators/NetworkName/staged_releases/<Data_Name>
#
# See README.md for the config format and a worked example per native space.
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <config.ini> [staging_root]" >&2
    exit 2
fi
CONFIG=$1
STAGING_ROOT_OVERRIDE=${2:-}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Parse + validate the config up front; parse_atlas_config.py exits non-zero with a
# clear message (unsupported space, Hard-only violation, count mismatches, etc.)
# before anything expensive runs. This also gives us CNC_SOURCE_FILES already
# naturally sorted (matching sort_files()'s natsorted() in the toolbox itself), so
# everything below passes explicit, pre-sorted file lists into MATLAB instead of
# letting dir() re-resolve and potentially re-sort them lexicographically.
#
# Checked explicitly (not just via `set -e`) because `eval "$(failing_cmd)"` does
# NOT abort under `set -e` when the failing command prints nothing to stdout --
# `eval ""` is itself a successful no-op, so a parse_atlas_config.py crash would
# otherwise silently fall through to confusing "unbound variable" errors below
# instead of the actual error already printed to stderr.
if ! CONFIG_SHELL=$(python3 "$SCRIPT_DIR/parse_atlas_config.py" --emit-shell "$CONFIG"); then
    echo "stage_atlas.sh: failed to parse/validate '$CONFIG' (see error above)" >&2
    exit 1
fi
eval "$CONFIG_SHELL"

STAGING_ROOT=${STAGING_ROOT_OVERRIDE:-/home/rkong/storage/ruby/data/Collaborators/NetworkName/staged_releases/${CNC_NAME}}
mkdir -p "$STAGING_ROOT"

echo "Staging atlas '$CNC_NAME' (category=$CNC_CATEGORY, native space=$CNC_SPACE, type=$CNC_TYPE_LOWER)"
echo "Staging root: $STAGING_ROOT"
echo "Source files (natural-sorted order): ${CNC_SOURCE_FILES[*]}"

# Build a MATLAB cell-array-of-strings literal, e.g. {'a','b'} or {}.
build_matlab_cell() {
    local arr=("$@")
    if [ ${#arr[@]} -eq 0 ]; then
        printf '{}'
        return
    fi
    local parts=() t
    for t in "${arr[@]}"; do
        parts+=("'${t}'")
    done
    local joined
    joined=$(IFS=,; echo "${parts[*]}")
    printf '{%s}' "$joined"
}
run_matlab() {
    local cmd=$1
    matlab -nodisplay -nosplash -r "try; ${cmd}; catch e; disp(getReport(e)); exit(1); end; exit(0);"
}

case "$CNC_SPACE" in
    FSLMNI2mm)
        MATLAB_FILES=$(build_matlab_cell "${CNC_SOURCE_FILES[@]}")
        run_matlab "addpath('${SCRIPT_DIR}'); proj_mni_fanout(${MATLAB_FILES}, '${CNC_CATEGORY}', '${CNC_TYPE_LOWER}', '${STAGING_ROOT}')"
        ;;

    fs_LR_32k|fsaverage6)
        # parse_atlas_config.py already enforced exactly one matched file for these spaces
        SOURCE_FILE="${CNC_SOURCE_FILES[0]}"
        run_matlab "addpath('${SCRIPT_DIR}'); proj_native_to_mni('${CNC_SPACE}', '${SOURCE_FILE}', '${CNC_NAME}', '${CNC_CATEGORY}', '${STAGING_ROOT}')"
        # proj_native_to_mni.m already wrote the other surface space + FSLMNI2mm;
        # only the remaining two volumetric spaces are needed here.
        bash "$SCRIPT_DIR/proj_mni_fanout_volumes.sh" \
            "${STAGING_ROOT}/atlases/FSLMNI2mm/${CNC_CATEGORY}/${CNC_NAME}.nii.gz" \
            "${CNC_CATEGORY}" "${CNC_TYPE_LOWER}" "${STAGING_ROOT}"
        ;;

    LairdColin2mm|ShenColin1mm)
        bash "$SCRIPT_DIR/proj_volume_native_to_mni.sh" \
            "${CNC_SPACE}" "${CNC_SOURCE_FILES_PATTERN}" "${CNC_CATEGORY}" "${CNC_TYPE_LOWER}" "${STAGING_ROOT}"
        # proj_volume_native_to_mni.sh preserves each file's basename, so the staged
        # FSLMNI2mm copies have the same names (and therefore the same natural-sorted
        # order) as CNC_SOURCE_FILES -- build that list explicitly rather than
        # re-globbing, to keep the order guaranteed identical to what was validated.
        STAGED_MNI_FILES=()
        for f in "${CNC_SOURCE_FILES[@]}"; do
            STAGED_MNI_FILES+=("${STAGING_ROOT}/atlases/FSLMNI2mm/${CNC_CATEGORY}/$(basename "$f")")
        done
        MATLAB_FILES=$(build_matlab_cell "${STAGED_MNI_FILES[@]}")
        # the generic fan-out also re-derives this atlas's OWN native space via a
        # volumetric round trip; we overwrite that below with the real original.
        run_matlab "addpath('${SCRIPT_DIR}'); proj_mni_fanout(${MATLAB_FILES}, '${CNC_CATEGORY}', '${CNC_TYPE_LOWER}', '${STAGING_ROOT}')"
        ;;

    *)
        echo "stage_atlas.sh: unsupported native space '$CNC_SPACE'" >&2
        exit 1
        ;;
esac

# Preserve true native-space fidelity. The generic MNI fan-out round-trips through
# FSLMNI2mm for every other space, including (for LairdColin2mm/ShenColin1mm-native)
# the atlas's OWN native space -- a round trip is lower fidelity than the atlas's
# real native file, which we already have, so overwrite it here.
if [ "$CNC_SPACE" != "FSLMNI2mm" ]; then
    NATIVE_OUT_DIR="${STAGING_ROOT}/atlases/${CNC_SPACE}/${CNC_CATEGORY}"
    mkdir -p "$NATIVE_OUT_DIR"
    if [ "$CNC_SPACE" == "fs_LR_32k" ] || [ "$CNC_SPACE" == "fsaverage6" ]; then
        # single file, renamed to match the atlas name (same convention
        # proj_native_to_mni.m already uses for the other surface space's output)
        EXT="${SOURCE_FILE##*.}"
        cp "$SOURCE_FILE" "${NATIVE_OUT_DIR}/${CNC_NAME}.${EXT}"
    else
        # possibly multiple files (e.g. a multi-map soft atlas, such as an ICA
        # decomposition with one file per component); preserve each original
        # basename, which must already contain the atlas name as a substring for
        # the toolbox's own name-matching to find it later.
        for f in "${CNC_SOURCE_FILES[@]}"; do
            cp "$f" "$NATIVE_OUT_DIR/"
        done
    fi
    echo "Restored authentic native-space file(s) into $NATIVE_OUT_DIR (overriding any round-tripped copy)"
fi

echo
echo "Staging complete: $STAGING_ROOT"
echo "Inspect the staged files under $STAGING_ROOT/atlases/, then run:"
echo "  python3 ${SCRIPT_DIR}/promote_atlas.py ${CONFIG} ${STAGING_ROOT}"
