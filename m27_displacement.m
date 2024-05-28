clear; clc; close all;

dir = '*';
dir_processing = sprintf('%s/03_Processing', dir);
dir_results = sprintf('%s/04_Results/27_displacement', dir);

subject = {'sub-01' 'sub-02' 'sub-03' 'sub-04' 'sub-05' 'sub-06' 'sub-07' 'sub-08' 'sub-09' 'sub-10' 'sub-11' 'sub-12' };
sequence = {'ovs_pf68' 'zoomit_pf68' 'zoomit_pf78' 'zoomit_pfno'};
runnumber = {'run1' 'run2'};
runfolder = {'02_run1' '03_run2'};

mean_displ_mcf = zeros(length(subject):length(sequence)*length(runnumber));
mean_displ_mcf_outl = zeros(length(subject):length(sequence)*length(runnumber));
mean_displ_sct = zeros(length(subject):length(sequence)*length(runnumber));
mean_displ_sct_outl = zeros(length(subject):length(sequence)*length(runnumber));

%%
for sub = 1:length(subject)
    for seq = 1:length(sequence)
        fprintf('%s %s \n', subject{sub}, sequence{seq})
        for run = 1:length(runnumber)
            txtfile = sprintf('%s/%s/func/derivatives/%s/%s/regressors/%s_%s_cr_mc_params_mcf.txt', dir_processing, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run});
            fileID = fopen(txtfile);
            params_mcf = fscanf(fileID,'%f%f%f%f%f%f', [6 Inf])';
            fclose(fileID);

            txtfile = sprintf('%s/%s/func/derivatives/%s/%s/regressors/%s_%s_cr_mc_params_sct.txt', dir_processing, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run});
            fileID = fopen(txtfile);
            fgetl(fileID);
            params_sct = fscanf(fileID,'%f%f', [2 Inf])';
            fclose(fileID);

            displacement_mcf = sqrt( (diff(params_mcf(:,4))).^2 + (diff(params_mcf(:,5))).^2 + (diff(params_mcf(:,6))).^2 );
            mean_displ_mcf(sub,(seq-1)*2+run) = mean(displacement_mcf);
            mean_displ_mcf_outl(sub,(seq-1)*2+run) = sum(displacement_mcf > 2);

            % displacement_sct = sqrt(params_sct(:,1).^2 + params_sct(:,2).^2);
            % mean_displ_sct(sub,(seq-1)*2+run) = mean(abs(diff(displacement_sct)));
            % mean_displ_sct_outl(sub,(seq-1)*2+run) = sum(displacement_sct > 2);
        end
    end
end

sum(mean_displ_mcf_outl, 2)
sum(mean_displ_mcf_outl, 'all')

%%
close all;
txtfile = sprintf('%s/sub-10/func/derivatives/zoomit_pfno/03_run2/regressors/zoomit_pfno_run2_cr_mc_params_mcf.txt', dir_processing);
fileID = fopen(txtfile);
params_mcf = fscanf(fileID,'%f%f%f%f%f%f', [6 Inf])';
fclose(fileID);

txtfile = sprintf('%s/sub-10/func/derivatives/zoomit_pfno/03_run2/regressors/zoomit_pfno_run2_cr_mc_params_sct.txt', dir_processing);
fileID = fopen(txtfile);
fgetl(fileID);
params_sct = fscanf(fileID,'%f%f', [2 Inf])';
fclose(fileID);

displacement_mcf = sqrt( (diff(params_mcf(:,4))).^2 + (diff(params_mcf(:,5))).^2 + (diff(params_mcf(:,6))).^2 );
displacement_sct = sqrt( (diff(params_sct(:,1))).^2 + (diff(params_sct(:,2))).^2 );

magnitude_mcf = sqrt( ((params_mcf(:,4))).^2 + ((params_mcf(:,5))).^2 + ((params_mcf(:,6))).^2 );


tiledlayout(2,2);
nexttile
plot(params_mcf(:,4:6))
ylim([-3 3])
legend('x', 'y', 'z')
title('mcf param')

nexttile
plot(displacement_mcf)
ylim([-3 3])
title('mcf abs disp')
hold on
plot(magnitude_mcf)
legend
hold off

clc
fprintf('mean displacement is %.2f \x00B1 %.2f \n', mean(displacement_mcf), std(displacement_mcf))
fprintf('%d framewise displacements are > %.2fmm \n', sum(displacement_mcf > (mean(displacement_mcf)+1.5*std(displacement_mcf))), mean(displacement_mcf)+1.5*std(displacement_mcf))
fprintf('%d framewise displacements are > 75%%+1.5IQR (%.2f)\n', sum(displacement_mcf > (prctile(displacement_mcf, 75) + 1.5*iqr(displacement_mcf))), prctile(displacement_mcf, 75) + 1.5*iqr(displacement_mcf))
fprintf('mean magnitude is %.2f \x00B1 %.2f \n', mean(magnitude_mcf), std(magnitude_mcf) )
fprintf('%d frames have a magnitude > %.2fmm \n', sum(magnitude_mcf > (mean(magnitude_mcf)+1.5*std(magnitude_mcf))), mean(magnitude_mcf)+1.5*std(magnitude_mcf))
fprintf('%d frames have a magnitude > 75%%+1.5IQR (%.2f)\n', sum(magnitude_mcf > (prctile(magnitude_mcf, 75) + 1.5*iqr(magnitude_mcf))), prctile(magnitude_mcf, 75) + 1.5*iqr(magnitude_mcf))

sum(( (magnitude_mcf > (mean(magnitude_mcf)+1*std(magnitude_mcf))) + ([0; displacement_mcf > (mean(displacement_mcf)+1*std(displacement_mcf))]) ) > 0)

nexttile
plot(params_sct)
legend('x', 'y')
ylim([-3 3])
title('sct param')

nexttile
plot(abs(diff(displacement_sct)))
ylim([-3 3])
title('sct abs disp')

%%
close all
params_x = squeeze(niftiread(sprintf('%s/%s/func/derivatives/%s/%s/regressors/%s_%s_cr_mc_params_x.nii', dir_processing, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run})));
params_y = squeeze(niftiread(sprintf('%s/%s/func/derivatives/%s/%s/regressors/%s_%s_cr_mc_params_y.nii', dir_processing, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run})));

framewise = sqrt((diff(params_x, 1, 2)).^2 + (diff(params_y, 1, 2)).^2);
plot(mean(framewise))
hold on

magnitude = sqrt(params_x.^2 + params_y.^2);
plot(mean(magnitude))
legend({'framewise' 'magnitude'})
