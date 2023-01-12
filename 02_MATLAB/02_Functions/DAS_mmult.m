function [img_arr, time] = DAS_mmult(dataset,indexmatrix,imagearea)

boolisgpu       = isgpuarray(dataset.rfdata) || isgpuarray(indexmatrix.M);
                if boolisgpu
aa              = gpuDevice;
                end
finalsize       = [length(imagearea.x_arr) length(imagearea.y_arr)];

tic
                %Sum and shape the final array
alpha           = dataset.rfdata;
img_arr         = indexmatrix.M*alpha(:);
                if boolisgpu
aa.             wait;
                end

img_arr         = reshape(img_arr,finalsize);                
time            = toc;

end