function img_arr = DAS_elementalLoop(dataset,delaymatrix)
%A function that uses a triple loop to perform DAS.

    %get sizes/preallocate
    [i_dim, j_dim, k_dim]  = size(delaymatrix);
    img_arr             = zeros(i_dim,j_dim);

    %delay and sum via triple loop to a 2D array (the most elemental way to do it)
    for ii = 1:i_dim
        for jj = 1:j_dim
            for kk = 1:k_dim
                img_arr(ii,jj) = img_arr(ii,jj) + dataset.rfdata(delaymatrix(ii,jj,kk),kk);
            end
        end
    end

end