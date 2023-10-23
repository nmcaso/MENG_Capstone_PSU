function img_arr = DAS_index(dataset, delaymatrix)
%img_arr = DAS_index(dataset, indexmatrix)
%
% Performs delay and sum by means of indexing timeseries data into a 3D
% delay matrix and then adding the results along the sensor dimension.
%
% Use the classes Dataset and IndexMatrix to generate the inputs.

    % this is the most vectorized version of the function on the CPU in MATLAB
    switch delaymatrix.type
        case "index"

            if dataset.n_frames > 1
                img_arr = zeros([size(delaymatrix.M, [1 2]) dataset.n_frames]);
                for ii = 1:size(img_arr, 3)
                    rff = dataset.rfdata(:,:,ii);
                    img_arr(:,:,ii) = sum(rff(delaymatrix.M), 3);
                end

            else
                img_arr     = sum(dataset.rfdata(delaymatrix.M),3);
            end
        
        
        otherwise
        error("Wrong DAS function for this type of index matrix!")
    end

end