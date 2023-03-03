function img_arr = DAS_mmult_gpu(dataset,delaymatrix,imagearea,cuda_device)
%this function is identical to DAS_mmult except for the title (so I can
%speed test in the comparison script easily.

finalsize       = [length(imagearea.x_arr) length(imagearea.y_arr)];        % Get the reshaping dimensions.
img_arr         = delaymatrix.M*dataset.rfdata(:);                          % Index by sparse matrix multiplication into a flat array
img_arr         = reshape(img_arr,finalsize);                               % Reshape the image to a 2D array
cuda_device.wait();
end