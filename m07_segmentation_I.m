
function m07_segmentation_I(file)
    V = spm_vol(file);
    Y = spm_read_vols(V);

    V.fname =sprintf('%s_pos.nii', file(1:length(file)-4));
    V.pinfo(2) = 0;
    Z = Y - 1*min(Y,[],'all') + 2;

    spm_write_vol(V,Z);

    disp('matlab done!')

end
