
function m07_segmentation_II_anat(file)
    [filepath, ~, ~] = fileparts(file);

    V = spm_vol(file);
    Y = spm_read_vols(V);

    W  = spm_vol(sprintf('%s/megre3d_cr.nii', filepath));
    W.fname =sprintf('%s_cor.nii', file(1:length(file)-4));
    W.pinfo(2) = 0;
    Z = flip(Y,1);

    spm_write_vol(W,Z);

    disp('matlab done!')

end
