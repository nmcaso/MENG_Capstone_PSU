%%
% DAS reconstruction template for circular array

load('Data_sets/leaf.mat') %data file

sound_speed = sos*1000; %convert to m/s

time_vec = rfdata'; %each transducer has a row

sens_vec = [x0;z0]*1e-3; %convert to meters

%imaging region in meters
x_vec = [-0.02:0.1e-3:0.02]';
y_vec = [-0.02:0.1e-3:0.02];

dt = 1/(fs*1e6);

pulse_index = 0; %setting t=0 value
sens_count = length(sens_vec);
t_max = size(time_vec,2); %maximum index number for the time delay matrix
total_pixels = length(x_vec)*length(y_vec);


"DAS reconstruction"
tic
DAS_image_circ = reshape(DAS_vectorized(x_vec,y_vec,time_vec,sens_vec,sound_speed,dt,pulse_index,sens_count,t_max,total_pixels),length(x_vec),length(y_vec));
toc
% The output from DAS_vectorized is a row vector equal to the total number
% of pixels. Reshape to convert the row vector to the original m x n vector

figure
mx = max(max(DAS_image_circ));
DAS_image_circ = DAS_image_circ/mx; %normalize
a = 3;
log_img = mx*(log10(1+a*DAS_image_circ./mx)./log10(1+a)); %increase contrast near zero

DAS_hilb_filt = MVHT(DAS_image_circ,1,'crop'); %filter for rotational interference

imagesc(y_vec,x_vec,flip(DAS_hilb_filt(end:-1:1,:),1))
title('DAS')
colormap(gray)
xlabel('m')
ylabel('m')
colorbar
set(gca, 'YDir','normal')