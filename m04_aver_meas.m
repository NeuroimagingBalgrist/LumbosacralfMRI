function m04_aver_meas (varargin)
% full argument list:
% file,N_echoes,N_meas,dummy_avgslr
% created by: G. David
% adapted by: C. Kündig
% changes:  -new argument jobfile_slr
%           -renaming of file/function from preproc_sc_megre_measuremnts to m04_aver_meas

jobfile_slr = varargin{5};

if nargin < 1
    % select first file
   [f,p] = uigetfile('*','Select the first NIfTI file of the measurement');
   p = p(1:end-1);
   ext = '.nii';
else
   file = varargin{1};
   [p,f,ext] = fileparts(file);
end

% extract the basename from filename
tmp = strsplit(deblank(f),'_');
basename = tmp{1};

fields = 2; % for test scans
% fields = 2:3; % for TASCI Balgrist
% fields = 2:4; % for the phase stabilization study
% fields = 2:3; % for TASCI Nottwil

for k = fields
    basename = [basename '_' tmp{k}];
end

% enter number of echoes
if nargin < 2
    N_echoes = -1;
    n = 0;
    while ~ismember(N_echoes,0:8)
        n = n+1;
        if n > 1
            'Incorrect input! The input has to be a non-negative integer in the range of 1-8.';
        end
        prompt = {
            ''
            ''
            'How many echoes are there? (for medic, enter 0)'
        };
        N_echoes = input([sprintf('%s\n\n',prompt{:}) 'Enter number of echoes: '],'s');
        N_echoes = str2num(N_echoes);
    end
else
    N_echoes = varargin{2};
end

% enter number of measurements
if nargin < 3
    meas_idx = 0;
    n = 0;
    while ~ismember(meas_idx,1:100)
        n = n+1;
        if n > 1
            'Incorrect input! The input has to be a non-negative integer.';
        end
        prompt = {
            ''
            ''
            'Which measurements are used for averaging?'
        };
        meas_idx = input([sprintf('%s\n\n',prompt{:}) 'Enter the measurements in square brackets: '],'s');
        meas_idx = str2num(meas_idx);
        N_meas = length(meas_idx);
    end
else
    meas_idx = varargin{3};
    N_meas = length(meas_idx);
end

% extract the indices of mesurements to be averaged
if N_meas == 1
    [a] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = a;
elseif N_meas == 2
    [a,b] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b];
elseif N_meas == 3
    [a,b,c] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c];
elseif N_meas == 4
    [a,b,c,d] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c,d];
elseif N_meas == 5
    [a,b,c,d,e] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c,d,e];
elseif N_meas == 6
    [a,b,c,d,e,f] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c,d,e,f];
elseif N_meas == 7
    [a,b,c,d,e,f,g] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c,d,e,f,g];
elseif N_meas == 8
    [a,b,c,d,e,f,g,h] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c,d,e,f,g,h];
elseif N_meas == 9
    [a,b,c,d,e,f,g,h,i] = feval(@(x) x{:}, num2cell(meas_idx));
    indices = [a,b,c,d,e,f,g,h,i];
end

% initialization
if N_echoes == 0 % medic
    V = spm_vol([p filesep basename '_meas-1' ext]);
    I_meas = zeros([V.dim, N_meas]);
end
if N_echoes > 0
    V = spm_vol([p filesep basename '_meas-1_echo-1' ext]);
    I_meas_c1 = zeros([V.dim, N_meas]);
end
if N_echoes > 1
    I_meas_c2 = zeros([V.dim, N_meas]);
    I_meas_c12 = zeros([V.dim, N_meas]);
end
if N_echoes > 2
    I_meas_c3 = zeros([V.dim, N_meas]);
    I_meas_c123 = zeros([V.dim, N_meas]);
end
if N_echoes > 3
    I_meas_c4 = zeros([V.dim, N_meas]);
    I_meas_c1234 = zeros([V.dim, N_meas]);
end
if N_echoes > 4
    I_meas_c5 = zeros([V.dim, N_meas]);
    I_meas_c12345 = zeros([V.dim, N_meas]);
end

% load in images to be averaged
for k = 1:N_meas
    if N_echoes == 0 % medic
        I_meas(:,:,:,k) = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) ext]));
    end
    if N_echoes > 0
        I_meas_c1(:,:,:,k) = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-1' ext]));
    end
    if N_echoes > 1
        I_meas_c2(:,:,:,k)  = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-2' ext]));
        I_meas_c12(:,:,:,k) = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-rms-12' ext]));
    end
    if N_echoes > 2
        I_meas_c3(:,:,:,k)   = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-3' ext]));
        I_meas_c123(:,:,:,k) = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-rms-123' ext]));
    end
    if N_echoes > 3
        I_meas_c4(:,:,:,k)    = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-4' ext]));
        I_meas_c1234(:,:,:,k) = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-rms-1234' ext]));
    end
    if N_echoes > 4
        I_meas_c5(:,:,:,k)     = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-5' ext]));
        I_meas_c12345(:,:,:,k) = spm_read_vols(spm_vol([p filesep basename '_meas-' num2str(indices(k)) '_echo-rms-12345' ext]));
    end
end

%% Average across measurements 1: ARITHMETIC AVERAGE
if N_echoes == 0 % medic
    I       = mean(I_meas,4);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '.nii']; spm_write_vol(V,I);
end
if N_echoes > 0
    I_c1    = mean(I_meas_c1,4);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-1.nii']; spm_write_vol(V,I_c1);
end
if N_echoes > 1
    I_c2    = mean(I_meas_c2,4);
    I_c12   = mean(I_meas_c12,4);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-2.nii']; spm_write_vol(V,I_c2);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-rms-12.nii']; spm_write_vol(V,I_c12);
end
if N_echoes > 2
    I_c3    = mean(I_meas_c3,4);
    I_c123  = mean(I_meas_c123,4);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-3.nii']; spm_write_vol(V,I_c3);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-rms-123.nii']; spm_write_vol(V,I_c123);
end
if N_echoes > 3
    I_c4    = mean(I_meas_c4,4);
    I_c1234 = mean(I_meas_c1234,4);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-4.nii']; spm_write_vol(V,I_c4);
    V.fname = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-rms-1234.nii']; spm_write_vol(V,I_c1234);
end
if N_echoes > 4
    I_c5     = mean(I_meas_c5,4);
    I_c12345 = mean(I_meas_c12345,4);
    V.fname  = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-5.nii']; spm_write_vol(V,I_c5);
    V.fname  = [p filesep basename '_meas-avg' num2str(N_meas) '_echo-rms-12345.nii']; spm_write_vol(V,I_c12345);
end

%% Average across measurements 1: MOTION CORRECTION (SERIAL LONGITUDINAL REGISTRATION)
if nargin < 4
    dummy_avgslr = -1;
    n = 0;
    while (isnumeric(dummy_avgslr) && ~ismember(dummy_avgslr,0:1))
        n = n+1;
        if n > 1
            'Incorrect input! The input has to be either 0 or 1.';
        end
        prompt = {
            ''
            ''
            'Do you want to perform motion correction? 0 for no, 1 for yes.'
        };
        dummy_avgslr = input([sprintf('%s\n\n',prompt{:}) 'Enter number (0 or 1): '],'s');
        dummy_avgslr = str2double(dummy_avgslr);
    end
else
    dummy_avgslr = varargin{4};
end

if dummy_avgslr
    if N_echoes == 0 % medic
        for j = 1:N_meas
            if N_meas > 0, filenames(j,:) = [p filesep basename '_meas-' num2str(indices(j)) ext]; end
        end

         % run SLR
        run_slr(filenames, N_meas, jobfile_slr);

        % rename output
        movefile([p filesep 'avg_' basename '_meas-1' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) ext]);
        V1 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) ext]);

        % remove empty slices
        remove_empty_slices(V1);
        clear filenames;
    end
    if N_echoes > 0
        for j = 1:N_meas
            if N_meas > 0, filenames(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-1' ext]; end
        end

        % run SLR
        run_slr(filenames, N_meas, jobfile_slr);

        % rename output
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-1' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-1' ext]);
        V1 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-1' ext]);

        % remove empty slices
        remove_empty_slices(V1);
        clear filenames;
    end
    if N_echoes > 1
        for j = 1:N_meas
            if N_meas > 0, filenames1(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-2' ext]; end
            if N_meas > 0, filenames2(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-rms-12' ext]; end
        end
        run_slr(filenames1, N_meas, jobfile_slr);
        run_slr(filenames2, N_meas, jobfile_slr);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-2' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-2' ext]);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-rms-12' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-12' ext]);
        V1 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-2' ext]);
        V2 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-12' ext]);
        remove_empty_slices(V1);
        remove_empty_slices(V2);
        clear filenames1 filenames2;
    end
    if N_echoes > 2
        for j = 1:N_meas
            if N_meas > 0, filenames1(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-3' ext]; end
            if N_meas > 0, filenames2(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-rms-123' ext]; end
        end
        run_slr(filenames1, N_meas, jobfile_slr);
        run_slr(filenames2, N_meas, jobfile_slr);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-3' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-3' ext]);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-rms-123' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-123' ext]);
        V1 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-3' ext]);
        V2 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-123' ext]);
        remove_empty_slices(V1);
        remove_empty_slices(V2);
        clear filenames1 filenames2;
     end
     if N_echoes > 3
        for j = 1:N_meas
            if N_meas > 0, filenames1(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-4' ext]; end
            if N_meas > 0, filenames2(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-rms-1234' ext]; end
        end
        run_slr(filenames1, N_meas, jobfile_slr);
        run_slr(filenames2, N_meas, jobfile_slr);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-4' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-4' ext]);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-rms-1234' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-1234' ext]);
        V1 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-4' ext]);
        V2 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-1234' ext]);
        remove_empty_slices(V1);
        remove_empty_slices(V2);
        clear filenames1 filenames2;
     end
    if N_echoes > 4
        for j = 1:N_meas
            if N_meas > 0, filenames1(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-5' ext]; end
            if N_meas > 0, filenames2(j,:) = [p filesep basename '_meas-' num2str(indices(j)) '_echo-rms-12345' ext]; end
        end
        run_slr(filenames1, N_meas, jobfile_slr);
        run_slr(filenames2, N_meas, jobfile_slr);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-5' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-5' ext]);
        movefile([p filesep 'avg_' basename '_meas-1' '_echo-rms-12345' ext],[p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-12345' ext]);
        V1 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-5' ext]);
        V2 = spm_vol([p filesep basename '_meas-avgslr' num2str(N_meas) '_echo-rms-12345' ext]);
        remove_empty_slices(V1);
        remove_empty_slices(V2);
        clear filenames1 filenames2;
    end
end

disp('Done!');

end

%% FUNCTIONS

% perform SLR
function run_slr(filenames, N_meas, jobfile)
    inputs = cell(2,1);
    inputs{1,1} = cellstr(filenames);
    inputs{2,1} = zeros(1,N_meas);

    % running job: serial longitudinal registration
    if ~isempty(filenames)
        spm('defaults','FMRI');
        spm_jobman('initcfg');
        spm_jobman('serial', jobfile, '', inputs{:});
    else
        disp('Subject did not run.');
    end
end


% remove empty slices
function remove_empty_slices(V1)
    I1 = spm_read_vols(V1);
    I1(I1>1e6) = 0;
    clear I1_new
    indices = zeros(V1.dim(3),1);
    for z = 1:V1.dim(3)
        if sum(sum(I1(:,:,z))) > 0, indices(z,1) = 1; end
    end
    I1_new = I1(:,:,indices==1);
    V1.dim = size(I1_new);
    V1.mat = V1.mat * spm_matrix([-1 -1 -1]*diag([1 1 1]))*spm_matrix([1,1,find(indices==1,1,'first')]);
    spm_write_vol(V1,I1_new);
end
