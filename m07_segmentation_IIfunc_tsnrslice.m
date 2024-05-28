
function m07_segmentation_IIfunc_tsnrslice(input, output)

    img = spm_vol(input);
    data = spm_read_vols(img);

    mean_slice = zeros(15,1);

    for i= 1:15
        mean_slice(i) = mean(nonzeros(data(:,:,i)));

    end
      
    writematrix(mean_slice,output)

    disp('matlab done!')

end
