function img_arr = DAS_elementalLoop(dataset,indexmatrix)
%A function that replicates the C elemental functionality of a triple loop
%to perform DAS.

    %get sizes/preallocate
    [xdim, ydim, zdim]  = size(indexmatrix);
    img_arr             = zeros(xdim,ydim);

    %delay and sum via triple loop to a 2D array (the most elemental way to do it)
    for ii = 1:xdim
        for jj = 1:ydim
            for kk = 1:zdim
                img_arr(ii,jj) = img_arr(ii,jj) + dataset.rfdata(indexmatrix(ii,jj,kk),kk);
            end
        end
    end

end