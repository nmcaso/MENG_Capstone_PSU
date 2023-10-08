function img_arr = DAS_half_vec_loop(dataset, delaymatrix)

    %preallocate zero array
    img_arr     = zeros(size(delaymatrix,[1 2]));
    
    %Add each page of the matrix
    for iter = 1:size(delaymatrix,3)

        rff         = squeeze(dataset.rfdata(:,iter));                      %extract a sensor's contrubution from rfdata
        img_arr     = img_arr + rff(delaymatrix(:,:,iter));                 %delay it and add it to the image

    end

end