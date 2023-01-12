function [img_arr, time] = DAS_parfor(dataset,indexmatrix)

tic

szs1to2     = size(indexmatrix.M);
img_arr     = zeros([szs1to2(1:2) 1]);

rf          = dataset.rfdata;
allind      = indexmatrix.M;

parfor iter = 1:szs1to2(3)
    
    rff         = squeeze(rf(:,iter));
    img_arr     = img_arr + rff(allind(:,:,iter));

end
time = toc;

end