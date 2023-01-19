function img_arr = DAS_half_vec_loop(dataset,indexmatrix)

    %preallocate zero array
    img_arr     = zeros(size(indexmatrix,[1 2]));

    %Add each page of the matrix
    for iter = 1:size(indexmatrix,3)

        rff         = squeeze(dataset.rfdata(:,iter));                      %extract a sensor's contrubution from rfdata
        img_arr     = img_arr + rff(indexmatrix(:,:,iter));                 %delay it and add it to the image
    
    end

end