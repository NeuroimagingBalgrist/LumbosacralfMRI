clear; clc; close all;

dir = '*';

subject = {'sub-01' 'sub-02' 'sub-03' 'sub-04' 'sub-05' 'sub-06' 'sub-07' 'sub-08' 'sub-09' 'sub-10' 'sub-11' 'sub-12'};
% subject = {'sub-01' 'sub-02' };
sequence = {'ovs_pf68' 'zoomit_pf68' 'zoomit_pf78' 'zoomit_pfno'};
runnumber = {'run1' 'run2'};
runfolder = {'02_run1' '03_run2'};

mean_epi_perc = zeros(4,12,2);
mean_epi_sm_perc = zeros(4,12,2);
mean_res_perc = zeros(4,12,2);


%%
for sub = 1:length(subject)
    for seq = 1:length(sequence)
        for run = 1:length(runnumber)
            epi = double(niftiread(sprintf('%s/%s/func/derivatives/%s/%s/%s_%s_cr_mc.nii', dir, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run})));
            epi_sm = double(niftiread(sprintf('%s/%s/func/derivatives/%s/%s/%s_%s_cr_mc_sm.nii', dir, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run})));
            res = double(niftiread(sprintf('%s/%s/func/derivatives/%s/%s/lvl1.feat/stats/res4d.nii.gz', dir, subject{sub}, sequence{seq}, runfolder{run})));
            func_mean = double(niftiread(sprintf('%s/%s/func/derivatives/%s/%s/lvl1.feat/mean_func.nii.gz', dir, subject{sub}, sequence{seq}, runfolder{run})));
            mask = niftiread(sprintf('%s/%s/func/derivatives/%s/%s/%s_%s_cr_mc_mrv.nii', dir, subject{sub}, sequence{seq}, runfolder{run}, sequence{seq}, runnumber{run}));

            epi_mean = mean(epi,4);
            epi_std = std(epi,0,4);
            epi_perc_mask = (epi_std./epi_mean.*100).*mask;
            mean_epi_perc(seq,sub,run) = mean(nonzeros(epi_perc_mask));

            epi_sm_mean = mean(epi_sm,4);
            epi_sm_std = std(epi_sm,0,4);
            epi_sm_perc_mask = (epi_sm_std./epi_sm_mean.*100).*mask;
            mean_epi_sm_perc(seq,sub,run) = mean(nonzeros(epi_sm_perc_mask));

            res_std = std(res,0,4);
            res_perc_mask = (res_std./func_mean.*100).*mask;
            mean_res_perc(seq,sub,run) = mean(nonzeros(res_perc_mask));

        end
    end
end

disp('std in percent for motioncorrected epi')
mean(mean_epi_perc,[2 3])
disp('std in percent for motioncorrected and smoothed epi')
mean(mean_epi_sm_perc,[2 3])
disp('std in percent for residual/mean_func')
mean(mean_res_perc,[2 3])
