function img_arr = DAS_mmult(dataset,delaymatrix,imagearea)
%image = DAS_mmult(dataset,indexmatrix,imagearea)
%
%Uses matrix multiplication to apply the appropriate time-delay to the
%timeseries' in dataset.rfdata using indexmatrix.M. This function also
%uses imagearea for the x and y image sizes (because indexmatrix.M is
%large and sparse, and not easy to tell the dimensions from).
%
%Use the Dataset, IndexMatrix, and ImageArea classes to generate the
%inputs.

if isa(delaymatrix, "DelayMatrix")
    img_arr         = delaymatrix.M*dataset.rfdata(:);                     % Index by sparse matrix multiplication into a flat array
else
    img_arr         = delaymatrix*dataset.rfdata(:);                       % same thing if the user put the M matrix itself into the delaymatrix input argument
end

img_arr         = reshape(img_arr,imagearea.image_size);                   % Reshape the image to a 2D array

end