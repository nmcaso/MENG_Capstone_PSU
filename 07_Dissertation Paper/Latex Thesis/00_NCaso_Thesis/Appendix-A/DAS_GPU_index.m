function img_arr = DAS_GPU_index(dataset,indexmatrix)
% img_arr = DAS_GPU_index(dataset, indexmatrix)
% 
% Creates a beamformed 2D Image by indexing on the GPU in 4 partitions

            %generate array indices in 4 quarters
q1          = 1:size(indexmatrix.M,3)/2;
q2          = size(indexmatrix.M,3)/4+1:size(indexmatrix.M,3)/2;
q3          = size(indexmatrix.M,3)/2+1:3*size(indexmatrix.M,3)/4;
q4          = 3*size(indexmatrix.M,3)/4+1:size(indexmatrix.M,3);

            %Sum and shape the final array
img_arr     = sum(dataset.rfdata(indexmatrix.M(:,:,q1)),3);
img_arr     = img_arr+sum(dataset.rfdata(indexmatrix.M(:,:,q2)),3);
img_arr     = img_arr+sum(dataset.rfdata(indexmatrix.M(:,:,q3)),3);
img_arr     = img_arr+sum(dataset.rfdata(indexmatrix.M(:,:,q4)),3);

end