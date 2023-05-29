function img_arr = DAS_index(dataset,delaymatrix)
%img_arr = DAS_index(dataset, indexmatrix)
%
% Performs delay and sum by means of indexing timeseries data into a 3D
% delay matrix and then adding the results along the sensor dimension.
%
% Use the classes Dataset and IndexMatrix to generate the inputs.

    % this is the most vectorized version of the function on the CPU in MATLAB
    switch delaymatrix.type
        case "Index2D"
        img_arr     = sum(dataset.rfdata(delaymatrix.M),3);
        otherwise
        error("Wrong DAS function for this type of index matrix!")
    end

end