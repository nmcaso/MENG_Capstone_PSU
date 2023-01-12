function [img_arr, time] = DAS_forloop(dataset,indexmatrix)

tic

szs1to2     = size(indexmatrix.M);
img_arr     = zeros([szs1to2(1:2) 1]);

for iter = 1:szs1to2(3)
    
    rff         = squeeze(dataset.rfdata(:,iter));
    img_arr     = img_arr + rff(indexmatrix.M(:,:,iter));

end
time = toc;

end