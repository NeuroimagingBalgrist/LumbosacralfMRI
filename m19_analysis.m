clearvars; clc; close all;

directory = '*';

%% Read Data

o68_path = sprintf('%s/grplvl_rando/00_results/ovs_pf68_tfce_corrp_tstat1_cr.nii', directory);
z68_path = sprintf('%s/grplvl_rando/00_results/zoomit_pf68_tfce_corrp_tstat1_cr.nii', directory);
z78_path = sprintf('%s/grplvl_rando/00_results/zoomit_pf78_tfce_corrp_tstat1_cr.nii', directory);
zNO_path = sprintf('%s/grplvl_rando/00_results/zoomit_pfno_tfce_corrp_tstat1_cr.nii', directory);
o68_head = spm_vol(o68_path); o68_image = spm_read_vols(o68_head);
z68_head = spm_vol(z68_path); z68_image = spm_read_vols(z68_head);
z78_head = spm_vol(z78_path); z78_image = spm_read_vols(z78_head);
zNO_head = spm_vol(zNO_path); zNO_image = spm_read_vols(zNO_head);
result_image = struct('o68_image', {o68_image}, 'z68_image', {z68_image}, 'z78_image', {z78_image}, 'zNO_image', {zNO_image});
result_name = fieldnames(result_image);


o68t_path = sprintf('%s/grplvl_rando/00_results/ovs_pf68_tstat1_cr.nii', directory);
z68t_path = sprintf('%s/grplvl_rando/00_results/zoomit_pf68_tstat1_cr.nii', directory);
z78t_path = sprintf('%s/grplvl_rando/00_results/zoomit_pf78_tstat1_cr.nii', directory);
zNOt_path = sprintf('%s/grplvl_rando/00_results/zoomit_pfno_tstat1_cr.nii', directory);
o68t_head = spm_vol(o68t_path); o68t_image = spm_read_vols(o68t_head);
z68t_head = spm_vol(z68t_path); z68t_image = spm_read_vols(z68t_head);
z78t_head = spm_vol(z78t_path); z78t_image = spm_read_vols(z78t_head);
zNOt_head = spm_vol(zNOt_path); zNOt_image = spm_read_vols(zNOt_head);
resultt_image = struct('o68t_image', {o68t_image}, 'z68t_image', {z68t_image}, 'z78t_image', {z78t_image}, 'zNOt_image', {zNOt_image});
resultt_name = fieldnames(resultt_image);

SC_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_cord.nii', directory);
SC_head = spm_vol(SC_path); SC_image = spm_read_vols(SC_head);

LD_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_LD_L3-S3.nii', directory);
RD_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_RD_L3-S3.nii', directory);
LV_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_LV_L3-S3.nii', directory);
RV_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_RV_L3-S3.nii', directory);
LVsMV_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_LV-MV_L3-S3.nii', directory);
RVsMV_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_RV-MV_L3-S3.nii', directory);
MV_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_region_MV_L3-S3.nii', directory);
LD_head = spm_vol(LD_path); LD_image = spm_read_vols(LD_head);
RD_head = spm_vol(RD_path); RD_image = spm_read_vols(RD_head);
LV_head = spm_vol(LV_path); LV_image = spm_read_vols(LV_head);
VR_head = spm_vol(RV_path); RV_image = spm_read_vols(VR_head);
LVsMV_head = spm_vol(LVsMV_path); LVsMV_image = spm_read_vols(LVsMV_head);
RVsMV_head = spm_vol(RVsMV_path); RVsMV_image = spm_read_vols(RVsMV_head);
MV_head = spm_vol(MV_path); MV_image = spm_read_vols(MV_head);
areaoi_image = struct('LD_image', {LD_image}, 'RD_image', {RD_image}, 'LV_image', {LV_image}, 'RV_image', {RV_image}, ...
    'LVsMV_image', {LVsMV_image}, 'RVsMV_image', {RVsMV_image},'MV_image', {MV_image});
areaoi_name = fieldnames(areaoi_image);


nL3_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_neuro_L3.nii', directory);
nL4_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_neuro_L4.nii', directory);
nL5_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_neuro_L5.nii', directory);
nS1_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_neuro_S1.nii', directory);
nS2_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_neuro_S2.nii', directory);
nS3_path = sprintf('%s/18_grplvl_mask/crop/PAM50_cr_neuro_S3.nii', directory);
nL3_head = spm_vol(nL3_path); nL3_image = spm_read_vols(nL3_head);
nL4_head = spm_vol(nL4_path); nL4_image = spm_read_vols(nL4_head);
nL5_head = spm_vol(nL5_path); nL5_image = spm_read_vols(nL5_head);
nS1_head = spm_vol(nS1_path); nS1_image = spm_read_vols(nS1_head);
nS2_head = spm_vol(nS2_path); nS2_image = spm_read_vols(nS2_head);
nS3_head = spm_vol(nS3_path); nS3_image = spm_read_vols(nS3_head);
nlevel_image = struct('nL3_image', {nL3_image}, 'nL4_image', {nL4_image}, 'nL5_image', {nL5_image}, 'nS1_image', {nS1_image}, ...
    'nS2_image', {nS2_image}, 'nS3_image', {nS3_image});
nlevel_name = fieldnames(nlevel_image);


%% Calculate Percent Activated per Level x Area

voxels_total = zeros(6, 7);
for k = 1:length(areaoi_name)
%     disp(areaoi_name{k})
    for l = 1:length(nlevel_name)
%         disp(nlevel_name{l})
        voxels_total(l, k) = sum(areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.99,"all");
    end
end

voxels_activated95 = zeros(6, 7, 4);
voxels_percent95 = zeros(6, 7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for k = 1:length(areaoi_name)
%         disp(areaoi_name{k})
        for l = 1:length(nlevel_name)
%             disp(nlevel_name{l})
            voxels_activated95(l, k, i) = sum(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.95,"all");
            voxels_percent95(l, k, i) = voxels_activated95(l, k, i)/voxels_total(l, k);
        end
    end
end

voxels_activated99 = zeros(6, 7, 4);
voxels_percent99 = zeros(6, 7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for k = 1:length(areaoi_name)
%         disp(areaoi_name{k})
        for l = 1:length(nlevel_name)
%             disp(nlevel_name{l})
            voxels_activated99(l, k, i) = sum(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.99,"all");
            voxels_percent99(l, k, i) = voxels_activated99(l, k, i)/voxels_total(l, k);
        end
    end
end

voxels_activated999 = zeros(6, 7, 4);
voxels_percent999 = zeros(6, 7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for k = 1:length(areaoi_name)
%         disp(areaoi_name{k})
        for l = 1:length(nlevel_name)
%             disp(nlevel_name{l})
            voxels_activated999(l, k, i) = sum(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.999,"all");
            voxels_percent999(l, k, i) = voxels_activated999(l, k, i)/voxels_total(l, k);
        end
    end
end

%% Calculate Percent Activated per Level or Area

nlevel_total = zeros(6, 4);
areaoi_total = zeros(7, 4);

nlevel_activated95  = zeros(6, 4);
nlevel_percent95  = zeros(6, 4);
areaoi_activated95  = zeros(7, 4);
areaoi_percent95  = zeros(7, 4);
nlevel_activated99  = zeros(6, 4);
nlevel_percent99  = zeros(6, 4);
areaoi_activated99  = zeros(7, 4);
areaoi_percent99  = zeros(7, 4);
nlevel_activated999  = zeros(6, 4);
nlevel_percent999  = zeros(6, 4);
areaoi_activated999  = zeros(7, 4);
areaoi_percent999  = zeros(7, 4);

for i = 1:length(result_name)
%     disp(result_name{i})
    for l = 1:length(nlevel_name)
%         disp(areaoi_name{l})
            nlevel_total(l, i) = sum(nlevel_image.(nlevel_name{l}),"all");
            nlevel_activated95(l, i) = sum(result_image.(result_name{i}).*nlevel_image.(nlevel_name{l})>0.95,"all");
            nlevel_percent95(l, i) = nlevel_activated95(l, i)/nlevel_total(l, i);
            nlevel_activated99(l, i) = sum(result_image.(result_name{i}).*nlevel_image.(nlevel_name{l})>0.99,"all");
            nlevel_percent99(l, i) = nlevel_activated99(l, i)/nlevel_total(l, i);
            nlevel_activated999(l, i) = sum(result_image.(result_name{i}).*nlevel_image.(nlevel_name{l})>0.999,"all");
            nlevel_percent999(l, i) = nlevel_activated999(l, i)/nlevel_total(l, i);
    end
%     Calculate percent activated per area from L3-S2 (S3 is empty).
    for k = 1:length(areaoi_name)
%          disp(nlevel_name{k})
            areaoi_total(k, i) = sum(areaoi_image.(areaoi_name{k})(:,:,33:121),"all");
            areaoi_activated95(k, i) = sum(result_image.(result_name{i})(:,:,33:121).*areaoi_image.(areaoi_name{k})(:,:,33:121)>0.95,"all");
            areaoi_percent95(k, i) = areaoi_activated95(k, i)/areaoi_total(k, i);
            areaoi_activated99(k, i) = sum(result_image.(result_name{i})(:,:,33:121).*areaoi_image.(areaoi_name{k})(:,:,33:121)>0.99,"all");
            areaoi_percent99(k, i) = areaoi_activated99(k, i)/areaoi_total(k, i);
            areaoi_activated999(k, i) = sum(result_image.(result_name{i})(:,:,33:121).*areaoi_image.(areaoi_name{k})(:,:,33:121)>0.999,"all");
            areaoi_percent999(k, i) = areaoi_activated999(k, i)/areaoi_total(k, i);
    end
end


%% Calculate Mean tScore per Level x Area in voxels aove threshold

tscore_mean95  = zeros(6, 7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for k = 1:length(areaoi_name)
%         disp(areaoi_name{k})
        for l = 1:length(nlevel_name)
%             disp(nlevel_name{l})
            tscore_mean95(l, k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.95)));
        end
    end
end

tscore_mean99  = zeros(6, 7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for k = 1:length(areaoi_name)
%         disp(areaoi_name{k})
        for l = 1:length(nlevel_name)
%             disp(nlevel_name{l})
            tscore_mean99(l, k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.99)));
        end
    end
end

tscore_mean999  = zeros(6, 7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for k = 1:length(areaoi_name)
%         disp(areaoi_name{k})
        for l = 1:length(nlevel_name)
%             disp(nlevel_name{l})
            tscore_mean999(l, k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k}).*nlevel_image.(nlevel_name{l})>0.999)));
        end
    end
end


%% Calculate Mean tScore per Level or Area in voxels above threshold

nlevel_tscore95  = zeros(6, 4);
areaoi_tscore95  = zeros(7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for l = 1:length(nlevel_name)
%         disp(areaoi_name{k})
            nlevel_tscore95(l, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*nlevel_image.(nlevel_name{l})>0.95)));
    end
    for k = 1:length(areaoi_name)
%             disp(nlevel_name{l})
            areaoi_tscore95(k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k})>0.95)));
    end
end

nlevel_tscore99  = zeros(6, 4);
areaoi_tscore99  = zeros(7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for l = 1:length(nlevel_name)
%         disp(areaoi_name{k})
            nlevel_tscore99(l, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*nlevel_image.(nlevel_name{l})>0.99)));
    end
    for k = 1:length(areaoi_name)
%             disp(nlevel_name{l})
            areaoi_tscore99(k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k})>0.99)));
    end
end

nlevel_tscore999  = zeros(6, 4);
areaoi_tscore999  = zeros(7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for l = 1:length(nlevel_name)
%         disp(areaoi_name{k})
            nlevel_tscore999(l, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*nlevel_image.(nlevel_name{l})>0.999)));
    end
    for k = 1:length(areaoi_name)
%             disp(nlevel_name{l})
            areaoi_tscore999(k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*(result_image.(result_name{i}).*areaoi_image.(areaoi_name{k})>0.999)));
    end
end


%% Calculate Mean tScore per Level or Area

nlevel_tscorefull  = zeros(6, 4);
areaoi_tscorefull  = zeros(7, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for l = 1:length(nlevel_name)
%         disp(areaoi_name{k})
            nlevel_tscorefull(l, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*nlevel_image.(nlevel_name{l})));
    end
    for k = 1:length(areaoi_name)
%             disp(nlevel_name{l})
            areaoi_tscorefull(k, i) = mean(nonzeros(resultt_image.(resultt_name{i}).*areaoi_image.(areaoi_name{k})));
    end
end


%% Figures 00_tscore.png

figure;
b = bar(nlevel_tscorefull(1:5,:));
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 5])
xticklabels({'L3','L4','L5','S1','S2'})
title({'Mean t-score per neurological level',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/00_tscore_neuro.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table(nlevel_tscorefull(1:5,:)','VariableNames', xticklabels)], sprintf('%s/19_analysis/00_tscore_neuro.csv', directory))


figure;
b = bar([areaoi_tscorefull(6,:); areaoi_tscorefull(2,:); areaoi_tscorefull(5,:); areaoi_tscorefull(1,:); areaoi_tscorefull(7,:)]);
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 5])
xticklabels({'RV','RD','LV','LD','MV'})
title({'Mean t-score per region',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/00_tscore_region.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table([areaoi_tscorefull(6,:); areaoi_tscorefull(2,:); areaoi_tscorefull(5,:); areaoi_tscorefull(1,:); areaoi_tscorefull(7,:)]', ...
    'VariableNames', xticklabels)], sprintf('%s/19_analysis/00_tscore_region.csv', directory))

%% Figures 99_percent.png

figure;
b = bar(nlevel_percent99(1:5,:));
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 1])
xticklabels({'L3','L4','L5','S1','S2'})
% set(gca,'YTick',[0,0.2,0.4,1])
title({'Ratio of p<0.01 per neurological level',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/99_percent_neuro.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table(nlevel_percent99(1:5,:)','VariableNames', xticklabels)], sprintf('%s/19_analysis/99_percent_neuro.csv', directory))


figure;
b = bar([areaoi_percent99(6,:); areaoi_percent99(2,:); areaoi_percent99(5,:); areaoi_percent99(1,:); areaoi_percent99(7,:)]);
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 1])
xticklabels({'RV','RD','LV','LD','MV'})
% set(gca,'YTick',[0,0.2,0.4,1])
title({'Ratio of p<0.01 per region',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/99_percent_region.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table([areaoi_percent99(6,:); areaoi_percent99(2,:); areaoi_percent99(5,:); areaoi_percent99(1,:); areaoi_percent99(7,:)]', ...
    'VariableNames', xticklabels)], sprintf('%s/19_analysis/99_percent_region.csv', directory))


%% Figures 95_percent.png

figure;
b = bar(nlevel_percent95(1:5,:));
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 1])
xticklabels({'L3','L4','L5','S1','S2'})
title({'Ratio of p<0.05 per neurological level',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/95_percent_neuro.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table(nlevel_percent95(1:5,:)','VariableNames', xticklabels)], sprintf('%s/19_analysis/95_percent_neuro.csv', directory))


figure;
b = bar([areaoi_percent95(6,:); areaoi_percent95(2,:); areaoi_percent95(5,:); areaoi_percent95(1,:); areaoi_percent95(7,:)]);
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 1])
xticklabels({'RV','RD','LV','LD','MV'})
title({'Ratio of p<0.05 per region',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/95_percent_region.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table([areaoi_percent95(6,:); areaoi_percent95(2,:); areaoi_percent95(5,:); areaoi_percent95(1,:); areaoi_percent95(7,:)]', ...
    'VariableNames', xticklabels)], sprintf('%s/19_analysis/95_percent_region.csv', directory))


%% Figures 999_percent.png

figure;
b = bar(nlevel_percent999(1:5,:));
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 1])
xticklabels({'L3','L4','L5','S1','S2'})
title({'Ratio of p<0.001 per neurological level',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/999_percent_neuro.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table(nlevel_percent999(1:5,:)','VariableNames', xticklabels)], sprintf('%s/19_analysis/999_percent_neuro.csv', directory))


figure;
b = bar([areaoi_percent999(6,:); areaoi_percent999(2,:); areaoi_percent999(5,:); areaoi_percent999(1,:); areaoi_percent999(7,:)]);
b(1).FaceColor = "#0072BD";
b(2).FaceColor = "#EDB120";
b(3).FaceColor = "#77AC30";
b(4).FaceColor = "#A2142F";
grid
ax = gca;
ax.FontSize = 25;
axis([0.5 5.5 0 1])
xticklabels({'RV','RD','LV','LD','MV'})
title({'Ratio of p<0.001 per region',' '}, 'FontSize', 21)
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/19_analysis/999_percent_region.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable([array2table(SEQ) array2table([areaoi_percent999(6,:); areaoi_percent999(2,:); areaoi_percent999(5,:); areaoi_percent999(1,:); areaoi_percent999(7,:)]', ...
    'VariableNames', xticklabels)], sprintf('%s/19_analysis/999_percent_region.csv', directory))


%% Calculate Histo tScore and Percent

tscore_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        tscore_histo(slice, i) = mean(nonzeros(resultt_image.(resultt_name{i})(:,:,slice).*SC_image(:,:,slice)));
    end
end


percent95_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        percent95_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*SC_image(:,:,slice)>0.95,'all')/sum(SC_image(:,:,slice),'all');
    end
end


percent99_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        percent99_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*SC_image(:,:,slice)>0.99,'all')/sum(SC_image(:,:,slice),'all');
    end
end


percent999_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        percent999_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*SC_image(:,:,slice)>0.999,'all')/sum(SC_image(:,:,slice),'all');
    end
end


%% Figure 99_hist_perc.png

slice = 33:1:121; %L3-S2 only

percent99_hist = figure;

plot(movmean(percent99_histo(33:121,1),10), slice, 'Color', "#0072BD", 'LineWidth',3.0)
hold on

plot(movmean(percent99_histo(33:121,2),10), slice, 'Color', "#EDB120", 'LineWidth',3.0)
hold on

plot(movmean(percent99_histo(33:121,3),10), slice, 'Color', "#77AC30", 'LineWidth',3.0)
hold on

plot(movmean(percent99_histo(33:121,4),10), slice, 'Color', "#A2142F", 'LineWidth',3.0)
axis([0 0.5 32 122])

ax = gca;
ax.FontSize = 25;
set(gca,'box','off')
ax.YAxis.Visible = 'off';
set(gcf, 'PaperUnits', 'centimeters');
x_width=15;y_width=40;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, sprintf('%s/19_analysis/99_percent_hist.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable(array2table(movmean(percent99_histo(33:121,:),10),'VariableNames', SEQ), sprintf('%s/19_analysis/99_percent_hist.csv', directory))

percent95_hist = figure;

plot(movmean(percent95_histo(33:121,1),10), slice, 'Color', "#0072BD", 'LineWidth',3.0)
hold on

plot(movmean(percent95_histo(33:121,2),10), slice, 'Color', "#EDB120", 'LineWidth',3.0)
hold on

plot(movmean(percent95_histo(33:121,3),10), slice, 'Color', "#77AC30", 'LineWidth',3.0)
hold on

plot(movmean(percent95_histo(33:121,4),10), slice, 'Color', "#A2142F", 'LineWidth',3.0)
axis([0 0.6 32 122])

ax = gca;
ax.FontSize = 25;
set(gca,'box','off')
ax.YAxis.Visible = 'off';
set(gcf, 'PaperUnits', 'centimeters');
x_width=15;y_width=40;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, sprintf('%s/19_analysis/95_percent_hist.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable(array2table(movmean(percent95_histo(33:121,:),10),'VariableNames', SEQ), sprintf('%s/19_analysis/95_percent_hist.csv', directory))

percent999_hist = figure;

plot(movmean(percent999_histo(33:121,1),10), slice, 'Color', "#0072BD", 'LineWidth',3.0)
hold on

plot(movmean(percent999_histo(33:121,2),10), slice, 'Color', "#EDB120", 'LineWidth',3.0)
hold on

plot(movmean(percent999_histo(33:121,3),10), slice, 'Color', "#77AC30", 'LineWidth',3.0)
hold on

plot(movmean(percent999_histo(33:121,4),10), slice, 'Color', "#A2142F", 'LineWidth',3.0)
axis([0 0.5 32 122])

ax = gca;
ax.FontSize = 25;
set(gca,'box','off')
ax.YAxis.Visible = 'off';
set(gcf, 'PaperUnits', 'centimeters');
x_width=15;y_width=40;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, sprintf('%s/19_analysis/999_percent_hist.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable(array2table(movmean(percent999_histo(33:121,:),10),'VariableNames', SEQ), sprintf('%s/19_analysis/999_percent_hist.csv', directory))


%% Figure 99_hist_perc_RV.png

percent99_RV_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        percent99_RV_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.RVsMV_image(:,:,slice)>0.99,'all')/sum(areaoi_image.RVsMV_image(:,:,slice),'all');
    end
end

slice = 33:1:121; %L3-S3 only

percent99_RV_hist = figure;

plot(movmean(percent99_RV_histo(33:121,1),10), slice, 'Color', "#0072BD", 'LineWidth',3.0)
hold on

plot(movmean(percent99_RV_histo(33:121,2),10), slice, 'Color', "#EDB120", 'LineWidth',3.0)
hold on

plot(movmean(percent99_RV_histo(33:121,3),10), slice, 'Color', "#77AC30", 'LineWidth',3.0)
hold on

plot(movmean(percent99_RV_histo(33:121,4),10), slice, 'Color', "#A2142F", 'LineWidth',3.0)
axis([0 1 32 122])

ax = gca;
ax.FontSize = 25;
set(gca,'box','off')
ax.YAxis.Visible = 'off';
set(gcf, 'PaperUnits', 'centimeters');
x_width=15;y_width=40;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, sprintf('%s/19_analysis/99_percent_hist_RV.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable(array2table(movmean(percent99_RV_histo(33:121,:),10),'VariableNames', SEQ), sprintf('%s/19_analysis/99_percent_hist_RV.csv', directory))

percent95_RV_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        percent95_RV_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.RVsMV_image(:,:,slice)>0.95,'all')/sum(areaoi_image.RVsMV_image(:,:,slice),'all');
    end
end


slice = 33:1:121; %L3-S3 only
percent95_RV_hist = figure;

plot(movmean(percent95_RV_histo(33:121,1),10), slice, 'Color', "#0072BD", 'LineWidth',3.0)
hold on

plot(movmean(percent95_RV_histo(33:121,2),10), slice, 'Color', "#EDB120", 'LineWidth',3.0)
hold on

plot(movmean(percent95_RV_histo(33:121,3),10), slice, 'Color', "#77AC30", 'LineWidth',3.0)
hold on

plot(movmean(percent95_RV_histo(33:121,4),10), slice, 'Color', "#A2142F", 'LineWidth',3.0)
axis([0 1 32 122])

ax = gca;
ax.FontSize = 25;
set(gca,'box','off')
ax.YAxis.Visible = 'off';
set(gcf, 'PaperUnits', 'centimeters');
x_width=15;y_width=40;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, sprintf('%s/19_analysis/95_percent_hist_RV.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable(array2table(movmean(percent95_RV_histo(33:121,:),10),'VariableNames', SEQ), sprintf('%s/19_analysis/95_percent_hist_RV.csv', directory))

percent999_RV_histo  = zeros(181, 4);
for i = 1:length(result_name)
%     disp(result_name{i})
    for slice = 1:181
        percent999_RV_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.RVsMV_image(:,:,slice)>0.999,'all')/sum(areaoi_image.RVsMV_image(:,:,slice),'all');
    end
end


slice = 33:1:121; %L3-S3 only
percent999_RV_hist = figure;

plot(movmean(percent999_RV_histo(33:121,1),10), slice, 'Color', "#0072BD", 'LineWidth',3.0)
hold on

plot(movmean(percent999_RV_histo(33:121,2),10), slice, 'Color', "#EDB120", 'LineWidth',3.0)
hold on

plot(movmean(percent999_RV_histo(33:121,3),10), slice, 'Color', "#77AC30", 'LineWidth',3.0)
hold on

plot(movmean(percent999_RV_histo(33:121,4),10), slice, 'Color', "#A2142F", 'LineWidth',3.0)
axis([0 1 32 122])

ax = gca;
ax.FontSize = 25;
set(gca,'box','off')
ax.YAxis.Visible = 'off';
set(gcf, 'PaperUnits', 'centimeters');
x_width=15;y_width=40;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, sprintf('%s/19_analysis/999_percent_hist_RV.png', directory))
% Save data as .csv
SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
writetable(array2table(movmean(percent999_RV_histo(33:121,:),10),'VariableNames', SEQ), sprintf('%s/19_analysis/999_percent_hist_RV.csv', directory))


%% Figure 99_hist_perc_LV.png

% percent99_LV_histo  = zeros(182, 4);
% for i = 1:length(result_name)
% %     disp(result_name{i})
%     for slice = 1:182
%         percent99_LV_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.LVsMV_image(:,:,slice)>0.99,'all')/sum(areaoi_image.LVsMV_image(:,:,slice),'all');
%     end
% end
%
% slice = 1:1:182;
%
% percent99_RV_hist = figure;
%
% plot(movmean(percent99_LV_histo(:,1),10), slice, 'Color', "#0072BD", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_LV_histo(:,2),10), slice, 'Color', "#EDB120", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_LV_histo(:,3),10), slice, 'Color', "#77AC30", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_LV_histo(:,4),10), slice, 'Color', "#A2142F", 'LineWidth',2.0)
% axis([0 1 0 114])
% grid
%
% ax = gca;
% ax.FontSize = 25;
% grid
% legend('OVS20', 'iFOV28', 'iFOV35', 'iFOV42', 'Location', 'NE')
%
% set(gcf, 'PaperUnits', 'centimeters');
% x_width=35;y_width=40;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf, sprintf('%s/19_analysis/percent99_LV_hist.png', directory))

%% Figure 99_hist_perc_RD.png

% percent99_RD_histo  = zeros(182, 4);
% for i = 1:length(result_name)
% %     disp(result_name{i})
%     for slice = 1:182
%         percent99_RD_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.RD_image(:,:,slice)>0.99,'all')/sum(areaoi_image.RD_image(:,:,slice),'all');
%     end
% end
%
% slice = 1:1:182;
%
% percent99_RV_hist = figure;
%
% plot(movmean(percent99_RD_histo(:,1),10), slice, 'Color', "#0072BD", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_RD_histo(:,2),10), slice, 'Color', "#EDB120", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_RD_histo(:,3),10), slice, 'Color', "#77AC30", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_RD_histo(:,4),10), slice, 'Color', "#A2142F", 'LineWidth',2.0)
% axis([0 1 0 114])
% grid
%
% ax = gca;
% ax.FontSize = 25;
% grid
% legend('OVS20', 'iFOV28', 'iFOV35', 'iFOV42', 'Location', 'NE')
%
% set(gcf, 'PaperUnits', 'centimeters');
% x_width=35;y_width=40;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf, sprintf('%s/19_analysis/percent99_RD_hist.png', directory))

%% Figure 99_hist_perc_LD.png

% percent99_LD_histo  = zeros(182, 4);
% for i = 1:length(result_name)
% %     disp(result_name{i})
%     for slice = 1:182
%         percent99_LD_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.LD_image(:,:,slice)>0.99,'all')/sum(areaoi_image.LD_image(:,:,slice),'all');
%     end
% end
%
% slice = 1:1:182;
%
% percent99_RV_hist = figure;
%
% plot(movmean(percent99_LD_histo(:,1),10), slice, 'Color', "#0072BD", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_LD_histo(:,2),10), slice, 'Color', "#EDB120", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_LD_histo(:,3),10), slice, 'Color', "#77AC30", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_LD_histo(:,4),10), slice, 'Color', "#A2142F", 'LineWidth',2.0)
% axis([0 1 0 114])
% grid
%
% ax = gca;
% ax.FontSize = 25;
% grid
% legend('OVS20', 'iFOV28', 'iFOV35', 'iFOV42', 'Location', 'NE')
%
% set(gcf, 'PaperUnits', 'centimeters');
% x_width=35;y_width=40;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf, sprintf('%s/19_analysis/percent99_LD_hist.png', directory))

%% Figure 99_hist_perc_MV.png

% percent99_MV_histo  = zeros(182, 4);
% for i = 1:length(result_name)
% %     disp(result_name{i})
%     for slice = 1:182
%         percent99_MV_histo(slice, i) = sum(result_image.(result_name{i})(:,:,slice).*areaoi_image.MV_image(:,:,slice)>0.99,'all')/sum(areaoi_image.MV_image(:,:,slice),'all');
%     end
% end
%
% slice = 1:1:182;
%
% percent99_MV_hist = figure;
%
% plot(movmean(percent99_MV_histo(:,1),10), slice, 'Color', "#0072BD", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_MV_histo(:,2),10), slice, 'Color', "#EDB120", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_MV_histo(:,3),10), slice, 'Color', "#77AC30", 'LineWidth',2.0)
% hold on
%
% plot(movmean(percent99_MV_histo(:,4),10), slice, 'Color', "#A2142F", 'LineWidth',2.0)
% axis([0 1 0 114])
% grid
%
% ax = gca;
% ax.FontSize = 25;
% grid
% legend('OVS20', 'iFOV28', 'iFOV35', 'iFOV42', 'Location', 'NE')
%
% set(gcf, 'PaperUnits', 'centimeters');
% x_width=35;y_width=40;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf, sprintf('%s/19_analysis/percent99_MV_hist.png', directory))

close all;
