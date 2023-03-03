function img_arr = DAS_original(dataset, img_area, sens_arr)
%The original function. Altered to use the class definitions generated
%later.
%  
% img_arr = DAS_original(dataset, img_area, sens_arr)
% 
% Performs Delay and Sum by generating a 2D index matrix for each sensor, mapping  
% the received timeseries signals for that sensor per the index matrix to a
% set of points, and then adding up the contributions for each sensor.

    % delay matrix for the first transducer (so we don't have to
    % preallocate)
    delay_mat = round(sqrt((img_area.x_arr-sens_arr.x0(1)).^2+(img_area.y_arr.'-sens_arr.z0(1)).^2)/dataset.c/dataset.dt);
    
    % correct bad indices
    delay_mat(delay_mat < 1) = 1;
    delay_mat(delay_mat > size(dataset.rfdata,1)) = size(dataset.rfdata,1);
    
    img_arr = dataset.rfdata(delay_mat,1); %Delay
    
    n_sens = size(dataset.rfdata,2);
    % repeat, adding the image for each transducer element
    for ii = 2:n_sens
    
        % generate index matrix
        delay_mat = round(sqrt((img_area.x_arr-sens_arr.x0(ii)).^2+(img_area.y_arr.'-sens_arr.z0(ii)).^2)/dataset.c/dataset.dt);
       
        %correct bad sensor indices
        delay_mat(delay_mat < 1) = 1;
        delay_mat(delay_mat > size(dataset.rfdata,1)) = size(dataset.rfdata,1);

        img_arr = img_arr + dataset.rfdata(delay_mat,ii); %Sum

    end

    img_arr = reshape(img_arr,[length(img_area.x_arr), length(img_area.y_arr)]);

end