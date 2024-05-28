
function m15_registration_PAM50(anat,txtfile)

    [filepath, filename, ~] = fileparts(anat);
    label = sprintf('%s/reg_temp/%s_label.nii', filepath, filename);

    fileID = fopen(txtfile);
    areadata = textscan(fileID,'%f%f','HeaderLines',3);
    fclose(fileID);

    % Get index of the maximum moving average of 3.
    means = movmean(areadata{1,2},3);
    [~,index] = max(means);

    [txtfilepath, ~]  = fileparts(txtfile);
    lsetxt = sprintf('%s/megre3d_cr_lse.txt', txtfilepath);
    writematrix(areadata{1,1}(index),lsetxt)

    % Create image matrix with 0 and label 59 at index slice for the LSE.
    image = zeros(97,97,20);
    image(49,49,areadata{1,1}(index)) = 59;

    % Get header from the megre3d and write with newly created matrix.
    A = spm_vol(anat);
    A.fname = label;

    spm_write_vol(A,image);

    disp('matlab done!')
end
