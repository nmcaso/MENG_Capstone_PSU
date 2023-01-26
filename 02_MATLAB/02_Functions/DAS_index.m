function img_arr = DAS_index(dataset,indexmatrix)
%img_arr = DAS_index(dataset, indexmatrix)
%
% Performs delay and sum by means of indexing timeseries data into a 3D
% delay matrix and then adding the results along the sensor dimension.
%
% Use the classes Dataset and IndexMatrix to generate the inputs.

            % this is the most vectorized version of the function on the CPU in MATLAB
switch indexmatrix.type
    case "Index2D"
    img_arr     = sum(dataset.rfdata(indexmatrix.M),3);
    case "Index3D"
    img_arr     = sum(dataset.rfdata(indexmatrix.M),4);
    otherwise
    error("Wrong DAS function for this type of index matrix!")
end