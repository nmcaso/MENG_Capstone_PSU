function img_arr = DAS_frame_rotate(dataset, imagearea, sensorarray, depth_correction, interpolate)

    lxar                = length(imagearea.x_arr);
    lx0                 = length(sensorarray.x0);
    lyar                = length(imagearea.y_arr);
    theta               = 180-180/pi*atan2(sensorarray.x0, sensorarray.z0);
    
    if lxar ~= lyar
        error("This method is only supported for square arrays currently")
    end

    frame_size = size(dataset.rfdata,1);
    Lx = abs(imagearea.x_arr(1)) + abs(imagearea.x_arr(end));
    Ly = abs(imagearea.y_arr(1)) + abs(imagearea.y_arr(end));
    dx = Lx/(lxar-1);
    dy = Ly/(lyar-1);

    expansion_factor = sqrt(2);

    inflated.x_arr = imagearea.x_arr(1)*expansion_factor:dx:imagearea.x_arr(end)*expansion_factor;
    inflated.y_arr = imagearea.y_arr(1)*expansion_factor:dy:imagearea.y_arr(end)*expansion_factor;
    
    if interpolate
    delay_matrix0 = depth_correction + sqrt((inflated.x_arr-sensorarray.x0(1)).^2+(inflated.y_arr.'-sensorarray.z0(1)).^2)/dataset.c/dataset.dt;
    else
    delay_matrix0 = uint16(depth_correction + sqrt((inflated.x_arr-sensorarray.x0(1)).^2+(inflated.y_arr.'-sensorarray.z0(1)).^2)/dataset.c/dataset.dt);
    end

    c_size = floor((size(delay_matrix0) - [lxar, lyar]) / 2);
    xind_vec = c_size(1):c_size(1) + lxar-1;
    yind_vec = c_size(2):c_size(2) + lyar-1;
    
    img_arr = dataset.rfdata(round(delay_matrix0(xind_vec,yind_vec)), 1);
    for ii = 2:lx0

        rotated_img = imrotate(delay_matrix0, theta(ii) ,'bilinear',"crop");
        delay_matrix = rotated_img(xind_vec, yind_vec);
        delay_matrix(delay_matrix < 1) = 1;
        delay_matrix(delay_matrix > frame_size) = frame_size;

        if interpolate
            delay_floor = floor(delay_matrix);
            delay_ceil = delay_floor + 1;
            delay_diffs = delay_matrix - delay_floor;

            img_arr = img_arr + dataset.rfdata(delay_floor, ii) ...
                + delay_diffs(:).*dataset.rfdata(delay_ceil,ii);
        else
            img_arr = img_arr + dataset.rfdata(delay_matrix,ii);
        end
        
    end

img_arr = reshape(img_arr,[lxar lyar]);
end