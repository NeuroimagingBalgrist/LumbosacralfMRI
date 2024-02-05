clear; clc;

dir = '*';

subject = {'sub-01' 'sub-02' 'sub-03' 'sub-04' 'sub-05' 'sub-06' 'sub-07' 'sub-08' 'sub-09' 'sub-10' 'sub-11' 'sub-12' };
sequence = {'ovs_pf68' 'zoomit_pf68' 'zoomit_pf78' 'zoomit_pfno'};

z_threshold = 2.3263;

percent_active = NaN(length(sequence), 2 * length(subject));
outliers = NaN(length(sequence), 2 * length(subject));


for seq = 1:length(sequence)

    a = 1;

    for sub = 1:length(subject)
        % Run1
        outliers_fID = fopen(sprintf('%s%s/func/derivatives/%s/02_run1/regressors/%s_run1_cr_mc_outliers_count.txt', dir, subject{sub}, sequence{seq}, sequence{seq}));
        outliers_formatSpec = '%d';
        outliers(seq,a) = fscanf(outliers_fID,outliers_formatSpec);

        mask_path = sprintf('%s%s/func/derivatives/%s/02_run1/effect_size/mask_rv.nii', dir, subject{sub}, sequence{seq});
        mask_head = spm_vol(mask_path);
        mask_img = spm_read_vols(mask_head);

        stats_path = sprintf('%s%s/func/derivatives/%s/02_run1/reg_anat/stats_old/zstat1.nii', dir, subject{sub}, sequence{seq});
        stats_head = spm_vol(stats_path);
        stats_img = spm_read_vols(stats_head);

        active_volume = stats_img > z_threshold;
        active_mask = active_volume .* mask_img;

        voxels_active = sum(active_mask, 'all');
        voxels_mask = sum(mask_img, 'all');

        percent_active(seq,a) = voxels_active / voxels_mask;
        a = a+1;

        % Run2
        outliers_fID = fopen(sprintf('%s%s/func/derivatives/%s/03_run2/regressors/%s_run2_cr_mc_outliers_count.txt', dir, subject{sub}, sequence{seq}, sequence{seq}));
        outliers_formatSpec = '%d';
        outliers(seq,a) = fscanf(outliers_fID,outliers_formatSpec);

        mask_path = sprintf('%s%s/func/derivatives/%s/03_run2/effect_size/mask_rv.nii', dir, subject{sub}, sequence{seq});
        mask_head = spm_vol(mask_path);
        mask_img = spm_read_vols(mask_head);

        stats_path = sprintf('%s%s/func/derivatives/%s/03_run2/reg_anat/stats_old/zstat1.nii', dir, subject{sub}, sequence{seq});
        stats_head = spm_vol(stats_path);
        stats_img = spm_read_vols(stats_head);


        active_volume = stats_img > z_threshold;
        active_mask = active_volume .* mask_img;

        voxels_active = sum(active_mask, 'all');
        voxels_mask = sum(mask_img, 'all');

        percent_active(seq,a) = voxels_active / voxels_mask;
        a = a+1;

    end
end

%% Create Figure
close all;
figure
tcl = tiledlayout(2,2);
title(tcl,sprintf('Outliers vs Ratio Voxels z > %.4f', z_threshold))

nexttile
scatter(outliers(1,:), percent_active(1,:))
title('OVS20')
grid on
axis([0 70 0 1])
xlabel('Outliers')
ylabel('Ratio active voxels')

nexttile
scatter(outliers(2,:), percent_active(2,:))
title('iFOV28')
grid on
axis([0 70 0 1])

nexttile
scatter(outliers(3,:), percent_active(3,:))
title('iFOV35')
grid on
axis([0 70 0 1])

nexttile
scatter(outliers(4,:), percent_active(4,:))
title('iFOV42')
grid on
axis([0 70 0 1])


set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'*/figures/outliersvsactivation.png')


%% Create Figure iFOV42
close all;
figure

nexttile
scatter(outliers(4,:), percent_active(4,:) ,30 ,'blue', 'filled')
grid on
axis([0 70 0 1])

title('iFOV42 Outliers vs Ratio of active voxels ', 'z > 2.33, 429 volumes/run')
xlabel('Outlier Volumes')
ylabel('Ratio')

ax = gca;
ax.FontSize = 20;
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'*/figures/outliersvsactivation_iFOV42.png')
