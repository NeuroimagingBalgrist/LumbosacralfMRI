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
        outliers_fID = fopen(sprintf('%s/03_Processing/%s/func/derivatives/%s/02_run1/regressors/%s_run1_cr_mc_outliers_count.txt', dir, subject{sub}, sequence{seq}, sequence{seq}));
        outliers_formatSpec = '%d';
        outliers(seq,a) = fscanf(outliers_fID,outliers_formatSpec);

        mask_path = sprintf('%s/03_Processing/%s/func/derivatives/%s/02_run1/effect_size/mask_rv.nii', dir, subject{sub}, sequence{seq});
        mask_head = spm_vol(mask_path);
        mask_img = spm_read_vols(mask_head);

        stats_path = sprintf('%s/03_Processing/%s/func/derivatives/%s/02_run1/reg_anat/stats_old/zstat1.nii', dir, subject{sub}, sequence{seq});
        stats_head = spm_vol(stats_path);
        stats_img = spm_read_vols(stats_head);

        active_volume = stats_img > z_threshold;
        active_mask = active_volume .* mask_img;

        voxels_active = sum(active_mask, 'all');
        voxels_mask = sum(mask_img, 'all');

        percent_active(seq,a) = voxels_active / voxels_mask;
        a = a+1;

        % Run2
        outliers_fID = fopen(sprintf('%s/03_Processing/%s/func/derivatives/%s/03_run2/regressors/%s_run2_cr_mc_outliers_count.txt', dir, subject{sub}, sequence{seq}, sequence{seq}));
        outliers_formatSpec = '%d';
        outliers(seq,a) = fscanf(outliers_fID,outliers_formatSpec);

        mask_path = sprintf('%s/03_Processing/%s/func/derivatives/%s/03_run2/effect_size/mask_rv.nii', dir, subject{sub}, sequence{seq});
        mask_head = spm_vol(mask_path);
        mask_img = spm_read_vols(mask_head);

        stats_path = sprintf('%s/03_Processing/%s/func/derivatives/%s/03_run2/reg_anat/stats_old/zstat1.nii', dir, subject{sub}, sequence{seq});
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
title(tcl, sprintf('Outliers vs Ratio Voxels z > %.2f', z_threshold))

nexttile
scatter(outliers(1,:)/4.11, percent_active(1,:))
title('OVS20')
grid on
axis([0 15 0 1])
xlabel('Outliers')
ylabel('Ratio active voxels')

nexttile
scatter(outliers(2,:)/5.05, percent_active(2,:))
title('iFOV28')
grid on
axis([0 15 0 1])

nexttile
scatter(outliers(3,:)/4.66, percent_active(3,:))
title('iFOV35')
grid on
axis([0 15 0 1])

nexttile
scatter(outliers(4,:)/4.29, percent_active(4,:))
title('iFOV42')
grid on
axis([0 15 0 1])


set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/24_outliers/outliersvsactivation.png', dir))


%% Create Figure iFOV42
close all;
figure

nexttile
scatter(outliers(4,:)/4.29, percent_active(4,:) ,30 ,'blue', 'filled')
grid on
axis([0 15 0 1])

title('iFOV42 Outliers vs Ratio of active voxels ', 'z > 2.33')
xlabel('Outliers [%]')
ylabel('Ratio')

ax = gca;
ax.FontSize = 20;
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/24_outliers/outliersvsactivation_iFOV42.png', dir))

%% Create Figure Sequence Comparison
close all
volumes = [411 505 466 429]';
outliers_perc = (outliers./volumes)*100;
outliers_perc_mean = mean(outliers_perc, 2);
outliers_perc_std = std(outliers_perc,0, 2);

plotcolors = hex2rgb(["#0072BD" "#EDB120" "#77AC30" "#A2142F"]);

figure;
boxplot(outliers_perc', 'Colors', plotcolors)
% lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% set(lines, 'Color', 'r');

set(findobj(gca,'type','line'),'linew', 1.5)
n = findobj(gcf,'tag','Outliers');
for j = 1:numel(n)
    n(j).MarkerEdgeColor = plotcolors(5-j,:);
end

axis([0.5 4.5 0 20])
ax = gca;
ax.FontSize = 20;
% xlabel('Sequence')
xticks([1,2,3,4])
xticklabels({"OVS20", "iFOV28", "iFOV35", "iFOV42"})
ylabel('Percent Outliers')
title('Percent Outliers per Sequence')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/24_outliers/sequencevsoutliers.png', dir))


%% All 8 runs ordered over time
order_seq = [4 1 2 3; 3 2 1 4; 2 3 4 1; 1 4 3 2; 1 4 2 3; 4 1 3 2; 2 3 1 4; 3 2 4 1; 3 1 2 4; 4 2 1 3; 2 4 3 1; 1 3 4 2;];
order_run = zeros(12,8);
order_run(:,1) = order_seq(:,1); order_run(:,3) = order_seq(:,2); order_run(:,5) = order_seq(:,3); order_run(:,7) = order_seq(:,4);
order_run(order_run==4) = 7; order_run(order_run==3) = 5; order_run(order_run==2) = 3;
order_run(:,2) = order_run(:,1)+1; order_run(:,4) = order_run(:,3)+1; order_run(:,6) = order_run(:,5)+1; order_run(:,8) = order_run(:,7)+1;

order_outliers = zeros(12,8);
index = 1;
for sub = 1:12
    order_outliers(sub,1) = outliers(order_seq(sub,:)==1,index);
    order_outliers(sub,2) = outliers(order_seq(sub,:)==1,index+1);
    order_outliers(sub,3) = outliers(order_seq(sub,:)==2,index);
    order_outliers(sub,4) = outliers(order_seq(sub,:)==2,index+1);
    order_outliers(sub,5) = outliers(order_seq(sub,:)==3,index);
    order_outliers(sub,6) = outliers(order_seq(sub,:)==3,index+1);
    order_outliers(sub,7) = outliers(order_seq(sub,:)==4,index);
    order_outliers(sub,8) = outliers(order_seq(sub,:)==4,index+1);

    index=index+2;
end

order_outliers_perc = zeros(12,8);
index = 1;
for sub = 1:12
    order_outliers_perc(sub,1) = outliers_perc(order_seq(sub,:)==1,index);
    order_outliers_perc(sub,2) = outliers_perc(order_seq(sub,:)==1,index+1);
    order_outliers_perc(sub,3) = outliers_perc(order_seq(sub,:)==2,index);
    order_outliers_perc(sub,4) = outliers_perc(order_seq(sub,:)==2,index+1);
    order_outliers_perc(sub,5) = outliers_perc(order_seq(sub,:)==3,index);
    order_outliers_perc(sub,6) = outliers_perc(order_seq(sub,:)==3,index+1);
    order_outliers_perc(sub,7) = outliers_perc(order_seq(sub,:)==4,index);
    order_outliers_perc(sub,8) = outliers_perc(order_seq(sub,:)==4,index+1);

    index=index+2;
end

close all
figure;
boxplot(order_outliers_perc, 'Colors', 'k')
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');

axis([0.5 8.5 0 20])
ax = gca;
ax.FontSize = 20;
xlabel('Run Number')
xticks([1,2,3,4,5,6,7,8])
xticklabels({"1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b"})
ylabel('Percent Outliers')
title('Percent Outliers over Time')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/24_outliers/runvsoutliers.png', dir))

%% Sequences ordered in time
close all

figure;
boxplot(reshape(order_outliers_perc,[24,4]), 'Colors', 'k')
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(findobj(gca,'type','line'),'linew', 1.5)

n = findobj(gcf,'tag','Outliers');
for j = 1:numel(n)
    n(j).MarkerEdgeColor = 'k';
end
axis([0.5 4.5 0 20])
ax = gca;
ax.FontSize = 20;
% xlabel('Sequence Number')
xticks([1,2,3,4])
xticklabels({"1", "2", "3", "4"})
ylabel('Percent Outliers')
title('Percent Outliers over Time')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/24_outliers/run_combvsoutliers.png', dir))

%% First vs second run for all sequences
order_outliers_run = [reshape(order_outliers_perc(:,[1 3 5 7]),[48,1]) reshape(order_outliers_perc(:,[2 4 6 8]),[48,1])];

close all
figure;
boxplot(order_outliers_run, 'Colors', 'k')
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');

axis([0.5 2.5 0 20])
ax = gca;
ax.FontSize = 20;
xlabel('Run Number per Sequence')
xticks([1,2])
xticklabels({"1", "2"})
ylabel('Percent Outliers')
title('Percent Outliers over Run-Rerun')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/24_outliers/run_rerunvsoutliers.png', dir))

%%
close all

outliers_reshaped = NaN(12:8);
for i = 1:12
    outliers_reshaped(i,:) = reshape(outliers(1:4,((i-1)*2+1):((i-1)*2+2))', [1 8]);
end

writematrix(outliers_reshaped,  sprintf('%s/04_Results/24_outliers/outliers_data.csv', dir))

for i = 1:4
    outliers_reshaped_perc(:,((i-1)*2+1):((i-1)*2+2)) = outliers_reshaped(:,((i-1)*2+1):((i-1)*2+2))/volumes(i)*100;
end

writematrix(outliers_reshaped_perc,  sprintf('%s/04_Results/24_outliers/outliers_perc_data.csv', dir))


outliers_perc_data = array2table(NaN(96, 6), "VariableNames", ["Outliers","Subject","Sequence","Order_S","Run","Order_R"]);
outliers_perc_data.Outliers = reshape(outliers_reshaped_perc', [96 1]);
outliers_perc_data.Subject = reshape(repmat(subject, [8 1]), [96 1]);
outliers_perc_data.Sequence = reshape(repmat({'OVS20' 'iFOV28' 'iFOV35' 'iFOV42'}, [2 12]), [96 1]);
outliers_perc_data.Order_S = reshape(repmat(reshape(order_seq', [48 1]), [1 2])', [96 1]);
outliers_perc_data.Run = repmat({'run1'; 'run2'}, [48 1]);
outliers_perc_data.Order_R = reshape(order_run', [96 1]);

writetable(outliers_perc_data, sprintf('%s/04_Results/24_outliers/outliers_perc_table.txt', dir))
