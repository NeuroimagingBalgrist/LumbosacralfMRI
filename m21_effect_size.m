clc; clearvars; close all;

%% Prepare variable and building blocks
dir = '*';
dir_processing = sprintf('%s/03_Processing', dir);
dir_results = sprintf('%s/04_Results/21_effect_size', dir);
subject = {'sub-01' 'sub-02' 'sub-03' 'sub-04' 'sub-05' 'sub-06' 'sub-07' 'sub-08' 'sub-09' 'sub-10' 'sub-11' 'sub-12' };
sequence = {'ovs_pf68' 'zoomit_pf68' 'zoomit_pf78' 'zoomit_pfno'};
run_dir = {'02_run1' '03_run2'};
TR_sequence = [1.460 1.190 1.290 1.400];

%% Read in image and recreate signal from GLM output

mean_rest = zeros(12,8);
mean_task = zeros(12,8);

for sub = 1:length(subject)
% for sub = 1:2
    disp(subject(sub))
    for seq = 1:length(sequence)
%     for seq = 1:1
        for run = 1:length(run_dir)
            TR = TR_sequence(seq);
            DIR_effectsize = sprintf('%s/%s/func/derivatives/%s/%s/effect_size', dir_processing, subject{sub}, sequence{seq}, run_dir{run});
            disp(DIR_effectsize)

            beta0_path  = sprintf('%s/mean_func.nii', DIR_effectsize);
            beta0_head  = spm_vol(beta0_path);
            beta0_image = spm_read_vols(beta0_head);

            beta1_path  = sprintf('%s/pe1.nii', DIR_effectsize);
            beta1_head  = spm_vol(beta1_path);
            beta1_image = spm_read_vols(beta1_head);

            mask_path   = sprintf('%s/mask_rv.nii', DIR_effectsize);
            mask_head   = spm_vol(mask_path);
            mask_image  = spm_read_vols(mask_head);

            res4d_path  = sprintf('%s/res4d.nii', DIR_effectsize);  % only for header and dimensions!
            res4d_head  = spm_vol(res4d_path);
            res4d_image = spm_read_vols(res4d_head);


            design_path = sprintf('%s/design.mat', DIR_effectsize); % paradigm
            design_table = readtable(design_path, 'FileType' , 'text','NumHeaderLines',5);
            design_matrix = table2array(design_table);


            expvar_image = zeros(size(res4d_image));    % image series containing the paradigm convoluted with the HRF
            expvar_array = squeeze(design_matrix(:,1));
            for vol = 1:size(res4d_image, 4)
                expvar_image(:,:,:,vol) = expvar_array(vol);
            end


            signal_image = zeros(size(res4d_image));    % reconstructed signal from GLM
            for vol = 1:size(signal_image,4)
                signal_image(:,:,:,vol) = (beta1_image(:,:,:).*expvar_image(:,:,:,vol) + beta0_image(:,:,:));
            end



            volumes = size(signal_image,4);
            mean_signal = zeros(1,volumes);

            for i = 1:volumes
                mean_signal(i) = mean(nonzeros(mask_image(:,:,:).*signal_image(:,:,:,i)));
            end

            % decide which volumes belong to rest and task while dropping
            % the first 3 volumes after a switch of condition (ramp of
            % BOLD)

            time = 0:TR:600;

            isodd = rem(floor((time)/15),2) == 1;
            istask = zeros(1,length(isodd));
            a = 2;
            istask(1) = isodd(1);
            while a <= length(isodd)
                if isodd(a) > isodd(a-1)
                    a = a+3;
                elseif isodd(a) == 1
                    istask(a) = 1;
                    a = a+1;
                else
                    a =a+1;
                end
            end

            iseven = ~isodd;
            isrest = zeros(1,length(iseven));
            a = 2;
            isrest(1) = iseven(1);
            while a <= length(iseven)
                if iseven(a) > iseven(a-1)
                    a = a+3;
                elseif iseven(a) == 1
                    isrest(a) = 1;
                    a = a+1;
                else
                    a =a+1;
                end
            end

            mean_task(sub,(seq-1)*2+run) = mean(nonzeros(istask.*mean_signal));
            mean_rest(sub,(seq-1)*2+run) = mean(nonzeros(isrest.*mean_signal));

        end
    end
end

%% Structure results in table, calculate effect size and output them

sequence_name = {'OVS20' 'iFOV28' 'iFOV35' 'iFOV42'};
results = table;
results.subject = subject';
for seq = 1:4
    for run  = 1:2
        results.(sprintf('%s-%s-rest', sequence_name{seq}, run_dir{run})) = mean_rest(:,(seq-1)*2+run);
        results.(sprintf('%s-%s-task', sequence_name{seq}, run_dir{run})) = mean_task(:,(seq-1)*2+run);
        results.(sprintf('%s-%s-percent', sequence_name{seq}, run_dir{run})) = 100*(mean_task(:,(seq-1)*2+run)./mean_rest(:,(seq-1)*2+run)-1);
    end
end

fileID = fopen(sprintf('%s/overview.txt', dir_results),'w');
fprintf(fileID, 'OVS20 mean = %0.3f, std = %0.3f\n', mean([results.("OVS20-02_run1-percent"); results.("OVS20-03_run2-percent")]), std([results.("OVS20-02_run1-percent"); results.("OVS20-03_run2-percent")]));
fprintf(fileID, 'iFOV28 mean = %0.3f, std = %0.3f\n', mean([results.("iFOV28-02_run1-percent"); results.("iFOV28-03_run2-percent")]), std([results.("iFOV28-02_run1-percent"); results.("iFOV28-03_run2-percent")]));
fprintf(fileID, 'iFOV35 mean = %0.3f, std = %0.3f\n', mean([results.("iFOV35-02_run1-percent"); results.("iFOV35-03_run2-percent")]), std([results.("iFOV35-02_run1-percent"); results.("iFOV35-03_run2-percent")]));
fprintf(fileID, 'iFOV42 mean = %0.3f, std = %0.3f\n', mean([results.("iFOV42-02_run1-percent"); results.("iFOV42-03_run2-percent")]), std([results.("iFOV42-02_run1-percent"); results.("iFOV42-03_run2-percent")]));
fclose(fileID);

writetable(results, sprintf('%s/results.csv', dir_results))
