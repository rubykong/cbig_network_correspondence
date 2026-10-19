function proj_mni_fanout(all_MNI_files, outsubdir, data_type, staging_root)
    % Generalizes the original proj_MNI_to_all.m: projects a native-FSLMNI2mm atlas
    % (one or more .nii.gz files) into fsaverage6, fs_LR_32k, and FSLMNI2mm, then
    % shells out to proj_mni_fanout_volumes.sh for LairdColin2mm/ShenColin1mm.
    %
    % Unlike the original, everything is written under
    % staging_root/atlases/<space>/outsubdir/ instead of directly into
    % cbig_network_correspondence_data, so the output can be inspected before
    % promote_atlas.py copies it into the real data repo.
    %
    % Data is carried through every projection as-is (no binarization here) --
    % a Soft/Metric atlas's continuous values are only ever thresholded at
    % runtime by the toolbox itself, via Data_Threshold in atlas_config (see
    % README.md). Binarizing before projection used to be done here so that
    % nearest-neighbor warping would have sharp boundaries, but Soft/Metric data
    % is warped with linear interpolation (see is_hard below), which blurs a
    % binarized map's edges straight back into floats -- so pre-projection
    % binarization never actually achieved sharp boundaries for the atlases it
    % was meant to help (e.g. ICA decompositions), it just discarded the
    % original continuous values for no benefit.
    %
    % all_MNI_files: cell array of explicit, already-sorted file paths (NOT a
    %              dir()-style glob -- the caller, stage_atlas.sh, resolves and
    %              naturally sorts these itself via parse_atlas_config.py/natsort,
    %              so that the order matches the order Network_Names_File's lines
    %              were validated against. Relying on dir()'s own ordering here
    %              previously risked a lexicographic-vs-natural-sort mismatch for
    %              unpadded multi-file atlases, e.g. an ICA decomposition named
    %              zstat1..zstat25.
    % outsubdir:   atlas category subdirectory name (same as Data_Category in the
    %              release config), e.g. 'Glasser'
    % data_type:   'hard', 'soft', or 'metric' (lowercase; controls interpolation)
    % staging_root: output root, e.g. Collaborators/NetworkName/staged_releases/<name>

    if isempty(all_MNI_files)
        error('proj_mni_fanout:no_files', 'all_MNI_files is empty');
    end

    fs6_outdir = fullfile(staging_root, 'atlases', 'fsaverage6', outsubdir);
    if ~exist(fs6_outdir, 'dir'); mkdir(fs6_outdir); end
    fs_LR_32k_outdir = fullfile(staging_root, 'atlases', 'fs_LR_32k', outsubdir);
    if ~exist(fs_LR_32k_outdir, 'dir'); mkdir(fs_LR_32k_outdir); end
    MNI_outdir = fullfile(staging_root, 'atlases', 'FSLMNI2mm', outsubdir);
    if ~exist(MNI_outdir, 'dir'); mkdir(MNI_outdir); end

    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'inflated', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'inflated', 'cortex');

    script_dir = fileparts(mfilename('fullpath'));
    volumes_script = fullfile(script_dir, 'proj_mni_fanout_volumes.sh');

    for i = 1:length(all_MNI_files)
        MNI_file = all_MNI_files{i};
        [~, MNI_file_name, ~] = fileparts(MNI_file);
        % fileparts only strips the last extension (.gz); drop the remaining .nii too
        MNI_file_name = MNI_file_name(1:end-4);

        tmp = MRIread(MNI_file);
        MRIwrite(tmp, fullfile(MNI_outdir, [MNI_file_name '.nii.gz']));

        % data_type here is 'hard'/'soft'/'metric' (the release config's
        % vocabulary, same as the toolbox's Data_Type). Soft and metric are both
        % continuous-valued and get linear interpolation; only hard parcellations
        % are nearest-neighbor. CBIG_project_fsaverage2fsLR below additionally
        % needs its own, different two-value vocabulary ('label'/'metric') --
        % passing 'hard'/'soft' straight through crashes it with "unknown type
        % of data", so translate explicitly rather than reusing data_type as-is.
        is_hard = strcmp(data_type, 'hard');
        if is_hard
            surf_interp = 'nearest';
            cbig_type_of_data = 'label';
        else
            surf_interp = 'linear';
            cbig_type_of_data = 'metric';
        end
        [lh_labels, rh_labels] = CBIG_ProjectMNI2fsaverage_Ants(tmp, 'fsaverage6', surf_interp);
        lh_labels(lh_avg_mesh.MARS_label == 1) = 0;
        rh_labels(rh_avg_mesh.MARS_label == 1) = 0;
        save(fullfile(fs6_outdir, [MNI_file_name '.mat']), 'lh_labels', 'rh_labels');

        % CBIG_project_fsaverage2fsLR refuses to run if folder_to_write already
        % exists (it rmdir's it on success); remove any stale leftover from an
        % earlier failed/interrupted run first.
        fsLR_tmp_dir = fullfile(fs_LR_32k_outdir, 'tmp');
        if exist(fsLR_tmp_dir, 'dir')
            rmdir(fsLR_tmp_dir, 's');
        end
        [lh_labels, rh_labels, ~, ~] = CBIG_project_fsaverage2fsLR( ...
            lh_labels, rh_labels, 'fsaverage6', cbig_type_of_data, fsLR_tmp_dir, '20170508');
        save(fullfile(fs_LR_32k_outdir, [MNI_file_name '.mat']), 'lh_labels', 'rh_labels');

        cmd = sprintf('bash %s %s %s %s %s', volumes_script, MNI_file, outsubdir, data_type, staging_root);
        status = system(cmd);
        if status ~= 0
            error('proj_mni_fanout:volumes_failed', ...
                'proj_mni_fanout_volumes.sh failed for %s (exit code %d)', MNI_file, status);
        end
    end
end
