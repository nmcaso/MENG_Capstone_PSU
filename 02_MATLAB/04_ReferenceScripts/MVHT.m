function img = MVHT(image,angle,style)
    angle_vec = 0:angle:360;
    filt_img = abs(hilbert(imrotate(image,angle_vec(1),style)));
    img = imrotate(filt_img,-angle_vec(1),style);
    for i=2:length(angle_vec)
        filt_img = abs(hilbert(imrotate(image,angle_vec(i),style)));
        img = img+imrotate(filt_img,-angle_vec(i),style);
    end
    img = img/length(angle_vec);
end

