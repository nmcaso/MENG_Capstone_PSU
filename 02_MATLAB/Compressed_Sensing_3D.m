            %% Image Compression Test
            clear; clc; close all;
            cd C:\Users\cason\OneDrive\Documents\PSU\Project\02_MATLAB
            for pth = strtrim(string(ls))'; if isfolder(pth); addpath(pth); end; end
           
n_avgs      = 50;

compression_ratio_time = 0.7;
compression_ratio_space = 0.5;

            %% Setup: Dataset, imagearea, img region

data        = Dataset("kwave_out3D.mat");
sens        = SensorArray("kwave_out3D.mat");
for fld     = ["x0" "z0"]
    sens.(fld) = cast(sens.(fld),"double");
end
sens.x0     = repmat(sens.x0(:), [1 length(sens.z0)]);
sens.z0     = repmat(sens.z0, [length(sens.z0) 1]);

            %% Do B-Mode FFT Method
padding_factor = 1;
rf_double   = [flip(data.rfdata,1); data.rfdata(2:end,:,:)];

k_t         = k_vector(repmat(rf_double, [padding_factor 1 1]) , data.c/data.fs/1e6);
k_t         = repmat(k_t, [1 sens.nSensors]);

dx          = mean(diff(sens.x0),'all');
k_x         = k_vector(sens.x0(:,1), dx);
k_x         = repmat(k_x', [padding_factor*length(rf_double) 1 sens.nSensors(2)]);

dz          = mean(diff(sens.z0),'all');
k_z(1,1,:)  = k_vector(sens.z0(1,:), dx);
k_z         = repmat(k_z, [padding_factor*length(rf_double) sens.nSensors(1) 1]);

k           = sqrt(k_t.^2 + k_x.^2 + k_z.^2);

omg_0       = data.c.*k_t;
omg_1       = data.c.*k;

scale       = data.c.^2 .* sqrt( (omg_0 ./ data.c).^2 - k_x.^2 - k_z.^2) ./ (2 .* omg_0);
scale       (omg_0 == 0 & k_x == 0 & k_z == 0) = data.c ./ 2;
nonhomogeneous = abs(omg_0) < data.c * hypot(k_x, k_z);

pad1         = [padding_factor*length(rf_double), sens.nSensors];
pad2         = [size(rf_double)];

tic
for ii = 1:n_avgs

    RF          = scale.*fftshift(fftn(ifftshift(rf_double), pad1));
    RF          (nonhomogeneous) = 0;

    RF_interp   = interp3(k_x, omg_0, k_z, RF, k_x, omg_1, k_z, "linear");
    RF_interp   (isnan(RF_interp)) = 0;
    p_recon     = real(fftshift(ifftn(ifftshift(RF_interp), pad2)));
end
time1       = toc;
time1       = time1/n_avgs;

p_recon     = flip(4/data.c.*permute(p_recon(data.frame_size:end,:,:), [2 3 1]),3);

            %% Compressive Sensing: B-Mode FFT Method
cs_inds.    time = false(2*data.frame_size, 1);
cs_inds.    space = false(sens.nSensors(1) ,1);

radius_t    = round((1-compression_ratio_time)/2 * size(rf_double,1));
center_t    = floor(size(rf_double,1)/2)+1;

radius_x    = round((1-compression_ratio_space)/2 * size(rf_double,2));
center_x    = floor(size(rf_double,2)/2)+1;

cs_inds     .time(center_t-radius_t:center_t+radius_t-1) = true;
cs_inds     .space(center_x - radius_x:center_x + radius_x - 1) = true;

k_x_sub     = k_x(cs_inds.time, cs_inds.space, cs_inds.space);
k_z_sub     = k_z(cs_inds.time, cs_inds.space, cs_inds.space);
omg_0_sub   = omg_0(cs_inds.time, cs_inds.space, cs_inds.space);
omg_1_sub   = omg_1(cs_inds.time, cs_inds.space, cs_inds.space);
nonhomogeneous_sub = nonhomogeneous(cs_inds.time, cs_inds.space, cs_inds.space);

M           = zeros(size(rf_double));

tic
for ii = 1:n_avgs 
    RF2         = scale.*fftshift(fftn(ifftshift(rf_double), pad1));
    RF2         = RF2(cs_inds.time, cs_inds.space, cs_inds.space);
    RF2         (nonhomogeneous_sub) = 0;

    RF_interp2  = interp3(k_x_sub, omg_0_sub, k_z_sub, RF2, k_x_sub, omg_1_sub, k_z_sub, "linear");
    RF_interp2  (isnan(RF_interp2)) = 0;
    M           (cs_inds.time, cs_inds.space, cs_inds.space) = RF_interp2;
    p_recon2     = real(fftshift(ifftn(ifftshift(M), pad2)));
end
time2       = toc;
time2       = time2/n_avgs;

p_recon2    = flip(4/data.c.*permute(p_recon2(data.frame_size:end,:,:), [2 3 1]), 3);

            %% Make Volume Plot of the reconstructed Image

max_intens  = max(abs(p_recon2),[],'all'); 
plot_cs     = p_recon2/max_intens;

max_intens2 = max(abs(p_recon), [], 'all');
plot_noncs  = p_recon/max_intens2;

num_intspts = 50;
inten_scale = max_intens/num_intspts:max_intens/num_intspts:max_intens;
alpha       = ((0:num_intspts-1)./num_intspts).^6;

queryPoints = linspace(min(inten_scale), max(inten_scale), 256);
alphamap    = interp1(inten_scale, alpha, queryPoints)';

fig3d1      = uifigure("Position",[100 100 500 1200]);
fig3d2      = uifigure("Position",[100+500+1 100 500 1200]);
fig3d3      = uifigure("Position",[100+500*2+2 100 500 1200]);

view        = viewer3d(BackgroundColor="white", BackgroundGradient=false, Lighting="off", Parent=fig3d1, CameraZoom=2, OrientationAxes="on");
cs_img      = volshow(plot_cs, Alphamap = alphamap, Parent=view, Colormap = pucolors.magma);

view2       = viewer3d(BackgroundColor="white", BackgroundGradient=false, Lighting="off", Parent=fig3d2, CameraZoom=2, OrientationAxes="on");
non_csimg   = volshow(plot_noncs, Alphamap = alphamap, Parent=view2, Colormap = pucolors.viridis);

view3       = viewer3d(BackgroundColor="white", BackgroundGradient=false, Lighting="off", Parent=fig3d3, CameraZoom=2, OrientationAxes="on");
diff_img    = volshow(abs(plot_noncs - plot_cs), Parent=view3, Colormap = pucolors.viridis);

exportapp(fig3d1, "05_DataOut\volplot3d_compressed.png");
exportapp(fig3d2, "05_DataOut\volplot3d_normal.png");
exportapp(fig3d3, "05_DataOut\volplot3d_diff.png");

            %% Back-Calculate the required pressure-field

M2          = zeros(size(rf_double));
M2          (cs_inds.time,cs_inds.space,cs_inds.space) = RF2;

p_backcalc  = real(fftshift(ifftn(ifftshift(M2))));
p_backcalc  = p_backcalc(data.frame_size+1:end,:);

speedup     =  ((time1 / time2) - 1) * 100;
error       = abs(p_recon - p_recon2)./abs(p_recon);
error       (isinf(error)) = 0;
pct_err     = geomean(error,'all')*100;

disp        ("Reconstruction Error: " + pct_err + "%")
disp        ("Non-compressed time: " + time1 + " ... Compressed time: " + time2)
disp        ("Speedup: " + speedup + "%")
disp        ("Final Compression Ratio: " +  100*(1-numel(RF2)/numel(RF)) + "%")
            %% Helper Functions

function vec_out = k_vector(vec_in, d_vecin)

    Nx = length(vec_in);

    if rem(Nx, 2) == 0
        % grid dimension has an even number of points
        nx = ((-Nx/2:Nx/2-1)/Nx).';
    else
        % grid dimension has an odd number of points
        nx = ((-(Nx-1)/2:(Nx-1)/2)/Nx).';
    end

    nx(floor(Nx/2) + 1) = 0;
    vec_out = (2*pi/d_vecin) .* nx;

end