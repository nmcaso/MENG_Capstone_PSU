function img_arr = DAS_GPU_index(dataset, delaymatrix)
% img_arr = DAS_GPU_index(dataset, indexmatrix)
% 
% Creates a beamformed 2D Image by indexing on the GPU in 4 partitions

            %generate array indices in 4 quarters
q1          = 1:size(delaymatrix,3)/2;
q2          = size(delaymatrix,3)/4+1:size(delaymatrix,3)/2;
q3          = size(delaymatrix,3)/2+1:3*size(delaymatrix,3)/4;
q4          = 3*size(delaymatrix,3)/4+1:size(delaymatrix,3);

            %Sum and shape the final array
img_arr     = sum(dataset.rfdata(delaymatrix(:,:,q1)),3);
img_arr     = img_arr+sum(dataset.rfdata(delaymatrix(:,:,q2)),3);
img_arr     = img_arr+sum(dataset.rfdata(delaymatrix(:,:,q3)),3);
img_arr     = img_arr+sum(dataset.rfdata(delaymatrix(:,:,q4)),3);

end