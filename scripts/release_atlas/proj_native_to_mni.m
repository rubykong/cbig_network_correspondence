function proj_native_to_mni(native_space, mat_file, atlas_name, outsubdir, staging_root)
    % Generalizes proj_Glasser_to_fs_MNI.m (native_space='fs_LR_32k') and
    % proj_Yeo_to_fslr_MNI.m (native_space='fsaverage6') into one function.
    %
    % Only supports Hard parcellations (a single lh_labels/rh_labels .mat) -- this
    % is the same restriction the two scripts it generalizes already had.
    %
    % Both native spaces need the atlas at full fsaverage7 resolution to call
    % CBIG_Projectfsaverage2MNI_Ants. fs_LR_32k gets this for free from
    % CBIG_project_fsLR2fsaverage; fsaverage6 needs an extra nearest-neighbor
    % upsample step (MARS_NNInterpolate_kdTree) since fsaverage6 vertices are the
    % first 40962 of fsaverage7's, not the reverse.
    %
    % native_space: 'fs_LR_32k' or 'fsaverage6'
    % mat_file:     path to the native .mat with lh_labels/rh_labels
    % atlas_name:   Data_Name, used for the output filename
    % outsubdir:    atlas category subdirectory (Data_Category), e.g. 'Glasser'
    % staging_root: output root, e.g. Collaborators/NetworkName/staged_releases/<name>
    %
    % Writes:
    %   staging_root/atlases/<other_surface_space>/outsubdir/<atlas_name>.mat
    %   staging_root/atlases/FSLMNI2mm/outsubdir/<atlas_name>.nii.gz
    % (the FSLMNI2mm file is then picked up by proj_mni_fanout.m/
    %  proj_mni_fanout_volumes.sh for the remaining LairdColin2mm/ShenColin1mm step)

    x = load(mat_file);
    % CBIG_project_fsLR2fsaverage/CBIG_project_fsaverage2fsLR below both require
    % their output folder to NOT already exist (they mkdir it themselves and
    % error "already exists" otherwise) -- so tmp_dir must only be removed if
    % stale from an earlier failed run, never pre-created.
    tmp_dir = fullfile(staging_root, 'tmp', atlas_name);
    if exist(tmp_dir, 'dir'); rmdir(tmp_dir, 's'); end

    switch native_space
        case 'fs_LR_32k'
            other_space = 'fsaverage6';
            other_outdir = fullfile(staging_root, 'atlases', other_space, outsubdir);
            if ~exist(other_outdir, 'dir'); mkdir(other_outdir); end

            [lh_FS7_data, rh_FS7_data] = CBIG_project_fsLR2fsaverage( ...
                x.lh_labels, x.rh_labels, 'fs_LR_32k', 'label', tmp_dir, '20170508');

            lh_labels = lh_FS7_data(1:40962); %#ok<NASGU>
            rh_labels = rh_FS7_data(1:40962); %#ok<NASGU>
            save(fullfile(other_outdir, [atlas_name '.mat']), 'lh_labels', 'rh_labels');

        case 'fsaverage6'
            other_space = 'fs_LR_32k';
            other_outdir = fullfile(staging_root, 'atlases', other_space, outsubdir);
            if ~exist(other_outdir, 'dir'); mkdir(other_outdir); end

            [lh_labels, rh_labels] = CBIG_project_fsaverage2fsLR( ...
                x.lh_labels, x.rh_labels, 'fsaverage6', 'label', tmp_dir, '20170508'); %#ok<ASGLU>
            save(fullfile(other_outdir, [atlas_name '.mat']), 'lh_labels', 'rh_labels');

            % upsample fsaverage6 -> fsaverage7 for the MNI projection below, since
            % fsaverage6 labels can't be exactly "densified" -- this is a nearest
            % neighbor approximation, same as the original proj_Yeo_to_fslr_MNI.m.
            lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
            rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex');
            lh_FS_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'sphere', 'cortex');
            rh_FS_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'sphere', 'cortex');

            lh_FS7_data = MARS_NNInterpolate_kdTree(lh_mesh7.vertices, lh_FS_mesh, x.lh_labels')';
            rh_FS7_data = MARS_NNInterpolate_kdTree(rh_mesh7.vertices, rh_FS_mesh, x.rh_labels')';

        otherwise
            error('proj_native_to_mni:unsupported_space', ...
                'native_space must be ''fs_LR_32k'' or ''fsaverage6'', got ''%s''', native_space);
    end

    mni_outdir = fullfile(staging_root, 'atlases', 'FSLMNI2mm', outsubdir);
    if ~exist(mni_outdir, 'dir'); mkdir(mni_outdir); end
    output = CBIG_Projectfsaverage2MNI_Ants(lh_FS7_data, rh_FS7_data, 'nearest');
    MRIwrite(output, fullfile(mni_outdir, [atlas_name '.nii.gz']));
end
