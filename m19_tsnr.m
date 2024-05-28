clearvars; clc; close all;

directory = '*';

subject = {'sub-01' 'sub-02' 'sub-03' 'sub-04' 'sub-05' 'sub-06' 'sub-07' 'sub-08' 'sub-09' 'sub-10' 'sub-11' 'sub-12' };
sequence = {'ovs_pf68' 'zoomit_pf68' 'zoomit_pf78' 'zoomit_pfno'};
subject_lse = zeros(1,12);
subject_tsnr = zeros(4,15,12);
subject_tsnr_csf = zeros(4,15,12);
subject_slice = zeros(15,12);


%% Get LSE slice for each subject
for subID = 1:length(subject)
    txtfile = sprintf('%s/03_Processing/%s/func/derivatives/segment_jim/megre3d_cr_lse.txt', directory, subject{subID});

    fileID = fopen(txtfile);
    subject_lse(subID) = fscanf(fileID,'%f')-5;
    fclose(fileID);

    subject_slice(:,subID) = 1-subject_lse(subID):1:15-subject_lse(subID);
end

%% Get tSNR for each subject for all sequences

for seqID = 1:length(sequence)
    for subID = 1:length(subject)
        txtfile = sprintf('%s/04_Results/sublvl_tsnr/%s_%s_rs_cr_mc_tsnr_slice_ero1.txt', directory, subject{subID}, sequence{seqID});

        fileID = fopen(txtfile);
        subject_tsnr(seqID,:,subID) = fscanf(fileID,'%f');
        fclose(fileID);

        txtfile = sprintf('%s/04_Results/sublvl_tsnr/%s_%s_rs_cr_mc_tsnr_csf_slice.txt', directory, subject{subID}, sequence{seqID});

        fileID = fopen(txtfile);
        subject_tsnr_csf(seqID,:,subID) = fscanf(fileID,'%f');
        fclose(fileID);
    end
end

%% Figures for each sequence

for seqID = 1:length(sequence)
    figure;
    plot(subject_slice,squeeze(subject_tsnr(seqID,:,:)))

    axis([-8 7 0 22])
    title({sprintf('%s',sequence{seqID}),' '}, 'FontSize', 21, 'Interpreter', 'none')
    set(gcf, 'PaperUnits', 'centimeters');
    x_width=20;y_width=15;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]);
    set (gca, 'xdir', 'reverse' )
    saveas(gcf, sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_%s.png', directory, sequence{seqID}))

end

for seqID = 1:length(sequence)
    figure;
    plot(subject_slice,squeeze(subject_tsnr_csf(seqID,:,:)))

    axis([-8 7 0 22])
    title({sprintf('%s',sequence{seqID}),' '}, 'FontSize', 21, 'Interpreter', 'none')
    set(gcf, 'PaperUnits', 'centimeters');
    x_width=20;y_width=15;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]);
    set (gca, 'xdir', 'reverse' )
    saveas(gcf, sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_csf_%s.png', directory, sequence{seqID}))

end

%% Figure mean and std

data = zeros(4,19,12);
for seqID = 1:length(sequence)
    for subID = 1:length(subject)
        offset = 13;
        slice = 1;
        for dataslice = 1:19
            if subject_slice(slice,subID) + offset == 1
                data(seqID,dataslice,subID) = subject_tsnr(seqID,slice,subID);
                if slice < 15
                    slice = slice + 1;
                end
            else
                data(seqID,dataslice,subID) = NaN;
            end
            offset = offset - 1;
        end
    end
end

data_mean = mean(data,3,"omitnan");
data_std = std(data,0,3,"omitnan");

figure;
% errorbar(-4.15:1:1.85,data_mean(1,9:15),data_std(1,9:15), Color = "#0072BD", LineWidth = 1.5)
plot(-4.15:1:1.85,data_mean(1,9:15), Color = "#0072BD", LineWidth = 1.5)
hold on
plot(-4.15:1:1.85,squeeze(data(1,9:15,:)), '+', Color = "#0072BD", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-4.05:1:1.95,data_mean(2,9:15),data_std(2,9:15), Color = "#EDB120", LineWidth = 1.5)
plot(-4.05:1:1.95,data_mean(2,9:15), Color = "#EDB120", LineWidth = 1.5)
hold on
plot(-4.05:1:1.95,squeeze(data(2,9:15,:)), '+', Color = "#EDB120", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-3.95:1:2.05,data_mean(3,9:15),data_std(3,9:15), Color = "#77AC30", LineWidth = 1.5)
plot(-3.95:1:2.05,data_mean(3,9:15), Color = "#77AC30", LineWidth = 1.5)
hold on
plot(-3.95:1:2.05,squeeze(data(3,9:15,:)), '+', Color = "#77AC30", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-3.85:1:2.15,data_mean(4,9:15),data_std(4,9:15), Color = "#A2142F", LineWidth = 1.5)
plot(-3.85:1:2.15,data_mean(4,9:15), Color = "#A2142F", LineWidth = 1.5)
hold on
plot(-3.85:1:2.15,squeeze(data(4,9:15,:)), '+', Color = "#A2142F", LineWidth = 1.0, MarkerSize = 3.0)

axis([-5 3 5 20])
ax = gca;
ax.FontSize = 20;
xlabel('Distance above LSE-slice [mm]')
xticks([-4,-3,-2,-1,0,1,2])
xticklabels({'-20','-15','-10','-5','0','5','10'})
ylabel('tSNR')
title({'Mean tSNR per sequences',' '})
% legend("OVS20", "iFOV28", "iFOV35", "iFOV42", 'Location', 'northeast')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'xdir', 'reverse' )
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_sc.png', directory))

%% Write data from figure to csv file

SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
column = {'-20' '-15' '-10' '-5' 'LSE' '5' '10'};
writetable([array2table(SEQ) array2table(data_mean(:,9:15),'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_sc_MEAN.csv', directory))
writetable([array2table(SEQ) array2table(data_std(:,9:15),'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_sc_STD.csv', directory))

close all;

%% Figure mean and std for csf

data_csf = zeros(4,19,12);
for seqID = 1:length(sequence)
    for subID = 1:length(subject)
        offset = 13;
        slice = 1;
        for dataslice = 1:19
            if subject_slice(slice,subID) + offset == 1
                data_csf(seqID,dataslice,subID) = subject_tsnr_csf(seqID,slice,subID);
                if slice < 15
                    slice = slice + 1;
                end
            else
                data_csf(seqID,dataslice,subID) = NaN;
            end
            offset = offset - 1;
        end
    end
end

data_csf_mean = mean(data_csf,3,"omitnan");
data_csf_std = std(data_csf,0,3,"omitnan");

figure;
% errorbar(-4.15:1:1.85,data_csf_mean(1,9:15),data_csf_std(1,9:15), Color = "#0072BD", LineWidth = 1.5)
plot(-4.15:1:1.85,data_csf_mean(1,9:15), Color = "#0072BD", LineWidth = 1.5)
hold on
plot(-4.15:1:1.85,squeeze(data_csf(1,9:15,:)), '+', Color = "#0072BD", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-4.05:1:1.95,data_csf_mean(2,9:15),data_csf_std(2,9:15), Color = "#EDB120", LineWidth = 1.5)
plot(-4.05:1:1.95,data_csf_mean(2,9:15), Color = "#EDB120", LineWidth = 1.5)
hold on
plot(-4.05:1:1.95,squeeze(data_csf(2,9:15,:)), '+', Color = "#EDB120", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-3.95:1:2.05,data_csf_mean(3,9:15),data_csf_std(3,9:15), Color = "#77AC30", LineWidth = 1.5)
plot(-3.95:1:2.05,data_csf_mean(3,9:15), Color = "#77AC30", LineWidth = 1.5)
hold on
plot(-3.95:1:2.05,squeeze(data_csf(3,9:15,:)), '+', Color = "#77AC30", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-3.85:1:2.15,data_csf_mean(4,9:15),data_csf_std(4,9:15), Color = "#A2142F", LineWidth = 1.5)
plot(-3.85:1:2.15,data_csf_mean(4,9:15), Color = "#A2142F", LineWidth = 1.5)
hold on
plot(-3.85:1:2.15,squeeze(data_csf(4,9:15,:)), '+', Color = "#A2142F", LineWidth = 1.0, MarkerSize = 3.0)

axis([-5 3 0 20])
ax = gca;
ax.FontSize = 20;
xlabel('Distance above LSE-slice [mm]')
xticks([-4,-3,-2,-1,0,1,2])
xticklabels({'-20','-15','-10','-5','0','5','10'})
ylabel('tSNR')
title({'Mean CSF tSNR per sequences',' '})
% legend("OVS20", "iFOV28", "iFOV35", "iFOV42", 'Location', 'northeast')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'xdir', 'reverse' )
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_csf.png', directory))

%% Write data from figure to csv file

SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
column = {'-20' '-15' '-10' '-5' 'LSE' '5' '10'};
writetable([array2table(SEQ) array2table(data_csf_mean(:,9:15),'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_csf_MEAN.csv', directory))
writetable([array2table(SEQ) array2table(data_csf_std(:,9:15),'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_csf_STD.csv', directory))

close all;

%% Ratio tSNR cord/csf

data_ratio = data./data_csf;

data_ratio_mean = mean(data_ratio,3,"omitnan");
data_ratio_std = std(data_ratio,0,3,"omitnan");

figure;
% errorbar(-4.15:1:1.85,data_ratio_mean(1,9:15),data_ratio_std(1,9:15), Color = "#0072BD", LineWidth = 1.5)
plot(-4.15:1:1.85,data_ratio_mean(1,9:15), Color = "#0072BD", LineWidth = 1.5)
hold on
plot(-4.15:1:1.85,squeeze(data_ratio(1,9:15,:)), '+', Color = "#0072BD", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-4.05:1:1.95,data_ratio_mean(2,9:15),data_ratio_std(2,9:15), Color = "#EDB120", LineWidth = 1.5)
plot(-4.05:1:1.95,data_ratio_mean(2,9:15), Color = "#EDB120", LineWidth = 1.5)
hold on
plot(-4.05:1:1.95,squeeze(data_ratio(2,9:15,:)), '+', Color = "#EDB120", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-3.95:1:2.05,data_ratio_mean(3,9:15),data_ratio_std(3,9:15), Color = "#77AC30", LineWidth = 1.5)
plot(-3.95:1:2.05,data_ratio_mean(3,9:15), Color = "#77AC30", LineWidth = 1.5)
hold on
plot(-3.95:1:2.05,squeeze(data_ratio(3,9:15,:)), '+', Color = "#77AC30", LineWidth = 1.0, MarkerSize = 3.0)
hold on

% errorbar(-3.85:1:2.15,data_ratio_mean(4,9:15),data_ratio_std(4,9:15), Color = "#A2142F", LineWidth = 1.5)
plot(-3.85:1:2.15,data_ratio_mean(4,9:15), Color = "#A2142F", LineWidth = 1.5)
hold on
plot(-3.85:1:2.15,squeeze(data_ratio(4,9:15,:)), '+', Color = "#A2142F", LineWidth = 1.0, MarkerSize = 3.0)

axis([-5 3 0 5])
ax = gca;
ax.FontSize = 20;
xlabel('Distance above LSE-slice [mm]')
xticks([-4,-3,-2,-1,0,1,2])
xticklabels({'-20','-15','-10','-5','0','5','10'})
ylabel('Ratio')
title({'Ratio Cord:CSF tSNR per sequences',' '})
% legend("OVS20", "iFOV28", "iFOV35", "iFOV42", 'Location', 'northeast')
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'xdir', 'reverse' )
set(gcf, 'PaperUnits', 'centimeters');
x_width=20;y_width=15;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf, sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_ratio.png', directory))

%% Write data from figure to csv file

SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
column = {'-20' '-15' '-10' '-5' 'LSE' '5' '10'};
writetable([array2table(SEQ) array2table(data_ratio_mean(:,9:15),'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_ratio_MEAN.csv', directory))
writetable([array2table(SEQ) array2table(data_ratio_std(:,9:15),'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_std_ratio_STD.csv', directory))

close all;

%% Write mean tSNR for cs, csf and ratio
data_results = NaN(4,6);

for i = 1:4
    data_results(i,1) = mean(data(i,9:15,:),'all',"omitnan");
    data_results(i,2) = std(data(i,9:15,:),0,'all',"omitnan");
    data_results(i,3) = mean(data_csf(i,9:15,:),'all',"omitnan");
    data_results(i,4) = std(data_csf(i,9:15,:),0,'all',"omitnan");
    data_results(i,5) = mean(data_ratio(i,9:15,:),'all',"omitnan");
    data_results(i,6) = std(data_ratio(i,9:15,:),0,'all',"omitnan");
end

SEQ = ["OVS20";"iFOV28";"iFOV35";"iFOV42"];
column = {'Cord' 'Cord_Std' 'CSF' 'CSF_Std' 'Cord/CSF' 'Cord/CSF_Std'};
writetable([array2table(SEQ) array2table(round(data_results,1), 'VariableNames', column)], sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_mean_results.csv', directory))

%% Repeated measurement ANOVA (done in R)
clc
ano_sub = 1:12;
ano_sub = repmat(ano_sub,28,1);
ano_sub = reshape(ano_sub,1,[])';

ano_seq = {'OVS20' 'iFOV28' 'iFOV35' 'iFOV42'}';
ano_seq = repmat(ano_seq,1,7);
ano_seq = reshape(ano_seq.',[],1);
ano_seq = repmat(ano_seq,12,1);

ano_slice = 1:7;
ano_slice = repmat(ano_slice,1,4);
ano_slice = repmat(ano_slice,1,12)';

ano_data = data(:,9:15,:);
ano_tsnr = reshape(permute(ano_data,[2 1 3]),[],1);

ano_table_slice = table(categorical(ano_sub), ano_seq, categorical(ano_slice), ano_tsnr, VariableNames={'Subject' 'Sequence' 'Slice' 'tSNR'});
writetable(ano_table_slice, sprintf('%s/04_Results/sublvl_tsnr/00_Analysis/tsnr_anova_data.txt', directory))

% aov_slice = anova(ano_table_slice, "tSNR ~ Sequence + Slice")
%
%
% ano_table_slice = table(categorical(ano_sub), ano_seq, categorical(ano_slice), ano_tsnr, VariableNames={'Subject' 'Sequence' 'Slice' 'tSNR'});
%
% % Assuming your dataset is stored in a table called 'data'
% % and the response variable is 'tSNR', and factors are 'Sequence', 'Slice', and 'Subject'
%
% % Create a repeated measures model
% rm = fitrm(ano_table_slice, 'tSNR ~ Sequence + Slice', 'WithinDesign', '');
%
% % Perform repeated measures ANOVA
% ranovaResults = ranova(rm, 'WithinModel', 'Slice*Subject');
%
% % Display the results
% disp(ranovaResults);
%
%
% % ano_sub = 1:12;
% % ano_sub = repmat(ano_sub,4,1);
% % ano_sub = reshape(ano_sub,1,[])';
% %
% % ano_seq = {'OVS20' 'iFOV28' 'iFOV35' 'iFOV42'}';
% % ano_seq = repmat(ano_seq,12,1);
% %
% % ano_data = mean(data(:,9:15,:),2);
% % ano_tsnr = reshape(permute(ano_data,[2 1 3]),[],1);
% %
% %
% % ano_table = table(categorical(ano_sub), ano_seq, ano_tsnr, VariableNames={'Subject' 'Sequence' 'tSNR'});
% %
% % aov = anova(ano_table, "tSNR ~ Sequence")
