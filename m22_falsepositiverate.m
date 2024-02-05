clear; clc;

dir = '*';

subject = {'sub-03n' 'sub-06n' 'sub-13n'};
sequence = {'ovs_pf68' 'zoomit_pf68' 'zoomit_pf78' 'zoomit_pfno'};

z_threshold = 2.3263;

fprate_volume = NaN(length(subject), 2*length(sequence));
fprate_mask = NaN(length(subject), 2*length(sequence));


for sub = 1:length(subject)
    a = 1;
    for seq = 1:length(sequence)
        mask_path = sprintf('%s%s/func/derivatives/%s/02_run1/%s_run1_cr_mc_mcsf_man.nii', dir, subject{sub}, sequence{seq}, sequence{seq});
        mask_head = spm_vol(mask_path);
        mask_img = spm_read_vols(mask_head);

        stats_path = sprintf('%s%s/func/derivatives/%s/02_run1/reg_anat/stats_old/zstat1.nii', dir, subject{sub}, sequence{seq});
        stats_head = spm_vol(stats_path);
        stats_img = spm_read_vols(stats_head);

        stats_pos = stats_img > z_threshold;

        pos_volume = sum(stats_pos,'all');
        vox_volume = stats_head.dim(1) * stats_head.dim(2) * stats_head.dim(3);
        fprate_volume(sub,a) = pos_volume/vox_volume;

        pos_mask = sum(stats_pos.*mask_img,'all');
        vox_mask = sum(mask_img,'all');
        fprate_mask(sub,a) = pos_mask/vox_mask;

        a = a+1;

        mask_path = sprintf('%s%s/func/derivatives/%s/03_run2/%s_run2_cr_mc_mcsf_man.nii', dir, subject{sub}, sequence{seq}, sequence{seq});
        mask_head = spm_vol(mask_path);
        mask_img = spm_read_vols(mask_head);

        stats_path = sprintf('%s%s/func/derivatives/%s/03_run2/reg_anat/stats_old/zstat1.nii', dir, subject{sub}, sequence{seq});
        stats_head = spm_vol(stats_path);
        stats_img = spm_read_vols(stats_head);

        stats_pos = stats_img > z_threshold;

        pos_volume = sum(stats_pos,'all');
        vox_volume = stats_head.dim(1) * stats_head.dim(2) * stats_head.dim(3);
        fprate_volume(sub,a) = pos_volume/vox_volume;

        pos_mask = sum(stats_pos.*mask_img,'all');
        vox_mask = sum(mask_img,'all');
        fprate_mask(sub,a) = pos_mask/vox_mask;

        a = a+1;
    end
end

fprintf('With z > %.4f the FP rate (volume) overall is %.4f +- %.4f \n', z_threshold, mean(fprate_volume,'all'), std(fprate_volume,0,'all'))
fprintf('With z > %.4f the FP rate (mask) overall is %.4f +- %.4f \n\n', z_threshold, mean(fprate_mask,'all'), std(fprate_mask,0,'all'))

fprintf('With z > %.4f the FP rate (volume) for %s is %.4f +- %.4f \n', z_threshold, sequence{1}, mean(fprate_volume(:,1:2),'all'), std(fprate_volume(:,1:2),0,'all'))
fprintf('With z > %.4f the FP rate (volume) for %s is %.4f +- %.4f \n', z_threshold, sequence{2}, mean(fprate_volume(:,3:4),'all'), std(fprate_volume(:,3:4),0,'all'))
fprintf('With z > %.4f the FP rate (volume) for %s is %.4f +- %.4f \n', z_threshold, sequence{3}, mean(fprate_volume(:,5:6),'all'), std(fprate_volume(:,5:6),0,'all'))
fprintf('With z > %.4f the FP rate (volume) for %s is %.4f +- %.4f \n', z_threshold, sequence{4}, mean(fprate_volume(:,7:8),'all'), std(fprate_volume(:,7:8),0,'all'))
