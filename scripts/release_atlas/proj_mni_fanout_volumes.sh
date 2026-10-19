#!/usr/bin/env bash
# Generalizes the original proj_MNI_to_diff_volume_space.sh: warps an
# already-staged FSLMNI2mm atlas file into LairdColin2mm and ShenColin1mm.
# Unlike the original, output goes under staging_root instead of directly into
# cbig_network_correspondence_data.
#
# Usage: proj_mni_fanout_volumes.sh <orig_file> <outsubdir> <data_type> <staging_root>
#   orig_file:    the native MNI file originally passed to proj_mni_fanout.m
#                 (only its basename is used, to locate the staged FSLMNI2mm copy)
#   outsubdir:    atlas category subdirectory (Data_Category), e.g. 'Glasser'
#   data_type:    'hard', 'soft', or 'metric' (lowercase; controls interpolation)
#   staging_root: same staging root passed to proj_mni_fanout.m
set -euo pipefail

orig_file=$1
outsubdir=$2
data_type=$3
staging_root=$4

# soft and metric are both continuous-valued and get linear interpolation; only
# hard parcellations are nearest-neighbor.
if [ "$data_type" == "hard" ]; then
    interp="NearestNeighbor"
else
    interp="Linear"
fi

filename=$(basename "$orig_file")
input="${staging_root}/atlases/FSLMNI2mm/${outsubdir}/${filename}"

# to LairdColin2mm
temp_2mm=/home/rkong/storage/ruby/data/Collaborators/NetworkName/Smith2009/maps_20ROIs/thresh_zstat1.nii.gz
warp_prefix=LairdColin2mm_to_FSLMNI2mm_fixed_ants
warpfolder=/home/rkong/storage/ruby/data/Collaborators/NetworkName/Registrations/LairdColin2mm_to_FSLMNI2mm
out_dir="${staging_root}/atlases/LairdColin2mm/${outsubdir}"
mkdir -p "$out_dir"
CBIG_antsApplyReg_vol2vol.sh -i "$input" -r "$temp_2mm" -w "$warp_prefix" -p "$filename" -d "$warpfolder" -o "$out_dir" -t "$interp" -s inverse

# to ShenColin1mm
temp_2mm=/home/rkong/storage/ruby/data/Collaborators/NetworkName/Shen/MNI_T1_1mm_stripped.nii.gz
warp_prefix=ShenColin1mm_to_FSLMNI2mm_fixed_ants
warpfolder=/home/rkong/storage/ruby/data/Collaborators/NetworkName/Registrations/ShenColin1mm_to_FSLMNI2mm
out_dir="${staging_root}/atlases/ShenColin1mm/${outsubdir}"
mkdir -p "$out_dir"
CBIG_antsApplyReg_vol2vol.sh -i "$input" -r "$temp_2mm" -w "$warp_prefix" -p "$filename" -d "$warpfolder" -o "$out_dir" -t "$interp" -s inverse
