function img_arr = DAS_vectorized(x_arr,y_arr,x0,z0,sens_arr,c,dt,pulse_ind,num_sens,t_max,pix_total)
    
    % initialize the delay matrix for the first transducer
    tic
    ind_mat = round(sqrt((x_arr-x0(20).').^2+(y_arr.'-z0(20)).^2)/c/dt+pulse_ind);
    toc
    
    %x_arr, y_arr, x0arr, y0arr, c, pulse_ind, c, sensnum
    
    %other:
    %an integer saying which
    
    %from sensorconfig:
    %x_arr, y_arr, pulse_ind
    
    %from dataset:
    %x0, z0, c, dt
    
    % check that the indeces are allowed
    if ind_mat<t_max & ind_mat>0
       img_arr = P_time(1,ind_mat);
       %Matlab lets you access a vector with a matrix of indeces. This
       %exploits that to generate the m x n image of pressures.
    else
       img_arr = zeros(1,pix_total);
    end
    % repeat, adding the image for each transducer element
    for i=2:num_sens
       tic
       ind_mat = round(sqrt((x_arr-sens_arr(1,i)).^2+(y_arr-sens_arr(2,i)).^2)/c/dt+pulse_ind);
       toc
       
       tic
       if ind_mat<t_max & ind_mat>0
          img_arr = img_arr+P_time(i,ind_mat);
       else
          img_arr = img_arr+zeros(1,pix_total);
       end
       toc
    end
end