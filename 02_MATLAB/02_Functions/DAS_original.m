function img_arr = DAS_original(dataset, img_area, sens_arr)
%The original function. Altered to use the class definitions generated
%later.
%  
% img_arr = DAS_original(img_area, sens_arr, dataset)
% 
% Performs Delay and Sum by generating a 2D index matrix for each sensor, mapping  
% the received timeseries signals for that sensor per the index matrix to a
% set of points, and then adding up the contributions for each sensor.

    % delay matrix for the first transducer (so we don't have to
    % preallocate)
    ind_mat = round(sqrt((img_area.x_arr-sens_arr.x0(1)).^2+(img_area.y_arr.'-sens_arr.z0(1)).^2)/dataset.c/dataset.dt);
    
    % correct bad indices
    ind_mat(ind_mat < 1) = 1;
    ind_mat(ind_mat > size(dataset.rfdata,1)) = size(dataset.rfdata,1);

    img_arr = dataset.rfdata(ind_mat,1); %Delay

    % repeat, adding the image for each transducer element
    for ii = 2:size(dataset.rfdata,2)
    
        % generate index matrix
        ind_mat = round(sqrt((img_area.x_arr-sens_arr.x0(ii)).^2+(img_area.y_arr.'-sens_arr.z0(ii)).^2)/dataset.c/dataset.dt);
       
        %correct bad sensor indices
        ind_mat(ind_mat < 1) = 1;
        ind_mat(ind_mat > size(dataset.rfdata,1)) = size(dataset.rfdata,1);

        img_arr = img_arr + dataset.rfdata(ind_mat,ii); %Sum

    end

    img_arr = reshape(img_arr,[length(img_area.x_arr), length(img_area.y_arr)]);

end