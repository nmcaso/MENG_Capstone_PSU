function DMAS_img = DMAS_gpu(x_arr,y_arr,P_time,sens_arr,c,dt,pulse_ind,num_sens, t_max,pix_total)
    ind_mat = round(sqrt(x_arr.^2+(y_arr-sens_arr(1)).^2)/c/dt+pulse_ind);
    img_arr = zeros(num_sens,pix_total);
 
    if ind_mat < t_max & ind_mat>0
       img_arr(1,:) = P_time(1,ind_mat);
    else
       img_arr(1,:) = zeros(1,pix_total);
    end

    for i=2:num_sens
       ind_mat = round(sqrt(x_arr.^2+(y_arr-sens_arr(i)).^2)/c/dt+pulse_ind);
       if ind_mat<t_max & ind_mat>0
          img_arr(i,:) = P_time(i,ind_mat);
       else
          img_arr(i,:) = zeros(1,pix_total);
       end
    end
    sign_mat = sign(img_arr);
    sqrt_mat = sqrt(abs(img_arr));
    DMAS_img = sum((sign_mat(2:end,:).*sqrt_mat(2:end,:)).*(sqrt_mat(1,:).*sign_mat(1,:)),1);
    for j=2:num_sens-1
        mult_img = sum((sign_mat(1:j-1,:).*sqrt_mat(1:j-1,:)).*(sqrt_mat(j,:).*sign_mat(j,:)),1);
        mult_img = mult_img+sum((sign_mat(j+1:end,:).*sqrt_mat(j+1:end,:).*(sqrt_mat(j,:).*sign_mat(j,:))),1);
        DMAS_img = DMAS_img+mult_img;
    end
    DMAS_img = DMAS_img+sum((sign_mat(1:-1,:).*sqrt_mat(1:-1,:).*(sqrt_mat(j,:).*sign_mat(j,:))),1);

end