function [img_arr, time] = DAS_index(dataset,indexmatrix)

tic
%Sum and shape the final array
img_arr         = mean(dataset.rfdata(indexmatrix.M),3);
time            = toc;

end