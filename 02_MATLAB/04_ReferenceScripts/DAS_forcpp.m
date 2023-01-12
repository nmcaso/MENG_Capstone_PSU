%% Full DAS including function

clear variables
load('Data_sets/leaf.mat') %data file

c = sos*1000; %convert to m/s

P_time = rfdata'; %each transducer has a row

sens_arr = [x0;z0]*1e-3; %convert to meters

%imaging region in meters
x_arr = [-0.02:0.1e-3:0.02].';
y_arr = [-0.02:0.1e-3:0.02];

dt = 1/(fs*1e6);

pulse_ind = 0; %setting t=0 value
num_sens = length(sens_arr);
t_max = size(P_time,2); %maximum index number for the time delay matrix
pix_total = length(x_arr)*length(y_arr);


fprintf('DAS Reconstruction: ');
tic
DAS_image_circ = reshape(DAS_vectorized(x_arr,y_arr,P_time,sens_arr,c,dt,pulse_ind,num_sens,t_max,pix_total),length(x_arr),length(y_arr));
toc
fprintf('\n');
% The output from DAS_vectorized is a row vector equal to the total number
% of pixels. Reshape to convert the row vector to the original m x n vector

imagesc(DAS_image_circ)

% function
function img_arr = DAS_vectorized(x_arr,y_arr,P_time,sens_arr,c,dt,pulse_ind,num_sens,t_max,pix_total)
    % initialize the delay matrix for the first transducer
    ind_mat = round(sqrt((x_arr-sens_arr(1,1)).^2+(y_arr-sens_arr(2,1)).^2)/c/dt+pulse_ind);
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
       ind_mat = round(sqrt((x_arr-sens_arr(1,i)).^2+(y_arr-sens_arr(2,i)).^2)/c/dt+pulse_ind);
       if ind_mat<t_max & ind_mat>0
          img_arr = img_arr+P_time(i,ind_mat);
       else
          img_arr = img_arr+zeros(1,pix_total);
       end
    end
end