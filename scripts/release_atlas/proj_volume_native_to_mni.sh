#!/usr/bin/env bash
# NEW (not a generalization of an existing script): forward-warps a
# LairdColin2mm-native or ShenColin1mm-native atlas into FSLMNI2mm, so it can then
# feed the existing MNI fan-out (proj_mni_fanout.m / proj_mni_fanout_volumes.sh)
# for everything else.
#
# This mirrors the forward-direction call pattern already used in production by
# apply_regAnts_Colin_to_MNI_Laird.sh (LairdColin2mm -> FSLMNI2mm, no "-s inverse"),
# just parametrized for either space and for an arbitrary set of input files instead
# of one hardcoded atlas. Nothing about the registration itself is new: same warp
# prefixes/folders/target template the rest of this pipeline already uses in the
# reverse ("-s inverse") direction for FSLMNI2mm -> LairdColin2mm/ShenColin1mm.
#
# Usage: proj_volume_native_to_mni.sh <native_space> <source_glob> <outsubdir> <data_type> <staging_root>
#   native_space: 'LairdColin2mm' or 'ShenColin1mm'
#   source_glob:  path or glob to the atlas's native-space .nii.gz file(s)
#   outsubdir:    atlas category subdirectory (Data_Category), e.g. 'Laird'
#   data_type:    'hard', 'soft', or 'metric' (lowercase; controls interpolation)
#   staging_root: output root, e.g. Collaborators/NetworkName/staged_releases/<name>
set -euo pipefail

native_space=$1
source_glob=$2
outsubdir=$3
data_type=$4
staging_root=$5

# soft and metric are both continuous-valued and get linear interpolation; only
# hard parcellations are nearest-neighbor.
if [ "$data_type" == "hard" ]; then
    interp="NearestNeighbor"
else
    interp="Linear"
fi

MNI_template=/home/rkong/storage/CBIG_private/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz

case "$native_space" in
    LairdColin2mm)
        warp_prefix=LairdColin2mm_to_FSLMNI2mm_fixed_ants
        warpfolder=/home/rkong/storage/ruby/data/Collaborators/NetworkName/Registrations/LairdColin2mm_to_FSLMNI2mm
        ;;
    ShenColin1mm)
        warp_prefix=ShenColin1mm_to_FSLMNI2mm_fixed_ants
        warpfolder=/home/rkong/storage/ruby/data/Collaborators/NetworkName/Registrations/ShenColin1mm_to_FSLMNI2mm
        ;;
    *)
        echo "proj_volume_native_to_mni.sh: unsupported native_space '$native_space' (must be LairdColin2mm or ShenColin1mm)" >&2
        exit 1
        ;;
esac

out_dir="${staging_root}/atlases/FSLMNI2mm/${outsubdir}"
mkdir -p "$out_dir"

shopt -s nullglob
files=( $source_glob )
shopt -u nullglob
if [ ${#files[@]} -eq 0 ]; then
    echo "proj_volume_native_to_mni.sh: no files matched '$source_glob'" >&2
    exit 1
fi

for input in "${files[@]}"; do
    filename=$(basename "$input")
    # forward direction (no "-s inverse"): native_space -> FSLMNI2mm
    CBIG_antsApplyReg_vol2vol.sh -i "$input" -r "$MNI_template" -w "$warp_prefix" \
        -p "$filename" -d "$warpfolder" -o "$out_dir" -t "$interp"
done
