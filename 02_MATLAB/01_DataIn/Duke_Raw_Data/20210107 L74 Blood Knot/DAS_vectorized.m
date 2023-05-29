function img_arr = robust_gpu(x_arr,y_arr,P_time,sens_arr,c,dt,pulse_ind,num_sens,t_max,pix_total)
    ind_mat = round(sqrt(x_arr.^2+(y_arr-sens_arr(1)).^2)/c/dt+pulse_ind);
    if ind_mat<t_max & ind_mat>0
       img_arr = P_time(1,ind_mat);
    else
       img_arr = zeros(1,pix_total);
    end
    for i=2:num_sens
       ind_mat = round(sqrt(x_arr.^2+(y_arr-sens_arr(i)).^2)/c/dt+pulse_ind);
       if ind_mat<t_max & ind_mat>0
          img_arr = img_arr+P_time(i,ind_mat);
       else
          img_arr = img_arr+zeros(1,pix_total);
       end
    end
end