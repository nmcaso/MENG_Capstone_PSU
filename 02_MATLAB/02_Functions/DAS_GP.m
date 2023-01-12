function [img_arr, time] = DAS_GP(dataset,indexmatrix)
aa = gpuDevice;
q1 = 1:size(indexmatrix.M,3)/2;
q2 = size(indexmatrix.M,3)/4+1:size(indexmatrix.M,3)/2;
q3 = size(indexmatrix.M,3)/2+1:3*size(indexmatrix.M,3)/4;
q4 = 3*size(indexmatrix.M,3)/4+1:size(indexmatrix.M,3);

tic
%Sum and shape the final array

img_arr         = mean(dataset.rfdata(indexmatrix.M(:,:,q1)),3);
img_arr         = img_arr+mean(dataset.rfdata(indexmatrix.M(:,:,q2)),3);
img_arr         = img_arr+mean(dataset.rfdata(indexmatrix.M(:,:,q3)),3);
img_arr         = img_arr+mean(dataset.rfdata(indexmatrix.M(:,:,q4)),3);
aa.wait();
time            = toc;

end