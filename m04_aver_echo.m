function m04_aver_echo (varargin)
% full argument list:
% file,N_echoes,N_meas
% created by: G. David
% adapted by: C. Kündig
% changes:  -renaming of file/function from preproc_sc_megre_echoes to m04_aver_echo

if nargin < 1
    % select first file
   [f,p] = uigetfile('*','Select the first NIfTI file of the measurement');
   p = p(1:end-1);
   ext = '.nii';
else
   file = varargin{1};
   [p,f,ext] = fileparts(file);
end

% extract basename from filename
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
    N_echoes = 0;
    n = 0;
    while ~ismember(N_echoes,1:8)
        n = n+1;
        if n > 1
            'Incorrect input! The input has to be a non-negative integer in the range of 1-8.';
        end
        prompt = {
            ''
            ''
            'How many echoes are there?'
        };
        N_echoes = input([sprintf('%s\n\n',prompt{:}) 'Enter number of echoes: '],'s');
        N_echoes = str2double(N_echoes);
    end
else
    N_echoes = varargin{2};
end

% enter number of measurements
if nargin < 3
    N_meas = 0;
    n = 0;
    while ~ismember(N_meas,1:8)
        n = n+1;
        if n > 1
            'Incorrect input! The input has to be a non-negative integer.';
        end
        prompt = {
            ''
            ''
            'How many measurements are there?'
        };
        N_meas = input([sprintf('%s\n\n',prompt{:}) 'Enter number of measurements: '],'s');
        N_meas = str2double(N_meas);
    end
else
    N_meas = varargin{3};
end

% combine echoes (simple average)
for k = 1:N_meas
    for i = 1:N_echoes
        filenames(i,:) = [p filesep basename '_meas-' num2str(k) '_echo-' num2str(i) ext];
    end
    Vall = spm_vol(filenames); V = Vall(1);
    I_meas = spm_read_vols(Vall);
    if N_echoes > 1
        I_meas_c12 = average_rms(I_meas(:,:,:,1:2));
        V.fname = [p filesep basename '_meas-' num2str(k) '_echo-rms-12' ext]; spm_write_vol(V,I_meas_c12);
    end
    if N_echoes > 2
        I_meas_c123 = average_rms(I_meas(:,:,:,1:3));
        V.fname = [p filesep basename '_meas-' num2str(k) '_echo-rms-123' ext]; spm_write_vol(V,I_meas_c123);
    end
    if N_echoes > 3
        I_meas_c1234 = average_rms(I_meas(:,:,:,1:4));
        V.fname = [p filesep basename '_meas-' num2str(k) '_echo-rms-1234' ext]; spm_write_vol(V,I_meas_c1234);
    end
    if N_echoes > 4
        I_meas_c12345 = average_rms(I_meas(:,:,:,1:5));
        V.fname = [p filesep basename '_meas-' num2str(k) '_echo-rms-12345' ext]; spm_write_vol(V,I_meas_c12345);
    end
    clear filenames I_meas;
end

disp('Done!');

end

% root mean square
function avg = average_rms (varargin)
    avg = 0;
    if nargin == 1
        dm = size(varargin{1});
        for i = 1:dm(4)
            avg = avg + varargin{1}(:,:,:,i).^2;
        end
        avg = sqrt(avg./dm(4));
    else
        for i = 1:nargin
            avg = avg + varargin{i}.^2;
        end
        avg = sqrt(avg./nargin);
    end
end

% arithmetic average
function avg = average_arithmetic (varargin)
    avg = 0;
    if nargin == 1
        dm = size(varargin{1});
        for i = 1:dm(4)
            avg = avg + varargin{i};
        end
        avg = sqrt(avg./dm(4));
    else
        for i = 1:nargin
            avg = avg + varargin{i};
        end
        avg = avg./nargin;
    end
end

% weighted average using squares as weights
function avg = average_wsq (varargin)
    avg = 0;
    sum = 0;
    if nargin == 1
        dm = size(varargin{1});
        for i = 1:dm(4)
            avg = avg + varargin{i}.^3;
            sum = sum + varargin{i}.^2;
        end
        avg = avg./sum;
    else
        for i = 1:nargin
            avg = avg + varargin{i}.^3;
            sum = sum + varargin{i}.^2;
        end
        avg = avg./sum;
    end
end
