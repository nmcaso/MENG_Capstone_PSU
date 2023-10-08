            %% Image Compression Test
            clear; clc; close all;
            for pth = strtrim(string(ls))'; if isfolder(pth); addpath(pth); end; end
           
n_avgs      = 1;
demo        = true;

compression_ratio_time = 0.8;
compression_ratio_space = 0.5;

            %% Demo!
if demo
    leaf        = imread("img.png").*8;
    kspace_leaf = dct2(leaf);
                tau = 2*pi;
    
    radius      = 440;
    x           = repmat(1:size(kspace_leaf,1),size(kspace_leaf,2),1);
    y           = repmat((1:size(kspace_leaf,2))',1,size(kspace_leaf,1));
    
    outside_radius = floor(hypot(x,y)) > radius;
    kspace_leaf(outside_radius) = 0;
    
    db_plot(kspace_leaf);
    title("Cosine Transform of Image (compressed)");
    leaf2       = abs(idct2(kspace_leaf));
    plot_two_leaves(leaf, leaf2, sum(kspace_leaf == 0,'all'));
end

            %% Setup: Dataset, imagearea, img region

data        = Dataset("kwave_out.mat");
sens        = SensorArray("kwave_out.mat");
for fld     = ["x0" "z0"]
    sens.(fld) = cast(sens.(fld),"double");
end

            %% Do B-Mode FFT Method
padding_factor = 1;
rf_double   = [flip(data.rfdata,1); data.rfdata];

k_t         = k_vector(repmat(rf_double,[padding_factor,1]) , data.c/data.fs);
k_t         = repmat(k_t, [1 sens.nSensors]);

dx          = mean(diff(sens.x0));
k_x         = k_vector(sens.x0, dx).';
k_x         = repmat(k_x, [padding_factor*length(rf_double) 1]);
k           = hypot(k_t, k_x);
omg_0       = data.c.*k_t;
omg_1       = data.c.*k;

pad1         = [padding_factor*length(rf_double), sens.nSensors];
pad2         = [2*data.frame_size, sens.nSensors];

tic
for ii = 1:n_avgs
    RF          = fftshift(fftn(ifftshift(rf_double), pad1));
    RF_interp   = interp2(k_x, omg_0, RF, k_x, omg_1, "cubic");
    RF_interp   (isnan(RF_interp)) = 0;
    p_recon     = real(fftshift(ifftn(ifftshift(RF_interp), pad2)));
end
time1       = toc;
time1       = time1/n_avgs;

p_recon     = p_recon(data.frame_size:end,:);

            %% Compressive Sensing: B-Mode FFT Method
cs_inds.    time = false(size(rf_double,1),1);
cs_inds.    space = false(size(rf_double,2),1);

radius_t    = round((1-compression_ratio_time)/2 * size(rf_double,1));
center_t    = floor(size(rf_double,1)/2)+1;

radius_x    = round((1-compression_ratio_space)/2 * size(rf_double,2));
center_x    = floor(size(rf_double,2)/2)+1;

cs_inds     .time(center_t-radius_t:center_t+radius_t-1) = true;
cs_inds     .space(center_x - radius_x:center_x + radius_x - 1) = true;

k_x_sub     = k_x(cs_inds.time,cs_inds.space);
omg_0_sub   = omg_0(cs_inds.time,cs_inds.space);
omg_1_sub   = omg_1(cs_inds.time,cs_inds.space);

M           = zeros(size(rf_double));

tic
for ii = 1:n_avgs 
    RF2         = fftshift(fftn(ifftshift(rf_double), pad1));
    RF2         = RF2(cs_inds.time,cs_inds.space);
    RF_interp2  = interp2(k_x_sub, omg_0_sub, RF2, k_x_sub, omg_1_sub, "linear");
    
    RF_interp2  (isnan(RF_interp2)) = 0;
    M           (cs_inds.time, cs_inds.space) = RF_interp2;
    p_recon2     = real(fftshift(ifftn(ifftshift(M), pad2)));
end
time2       = toc;
time2       = time2/n_avgs;

p_recon2    = p_recon2(data.frame_size:end,:);

figure      ("Position",[400 100 1600 900])
tiledlayout (1,3); 
nexttile;
imagesc     (p_recon2); colorbar;
title       (numel(p_recon2)/time2*1e-6 + " MegaPixels per sec at " + 100*(1-numel(RF2)/numel(RF))  + "% Compression");
subtitle    ("Speedup: " + time1/time2 + " \times");

nexttile;
imagesc     (20*log10(abs(p_recon2 - p_recon))./max(abs(p_recon2 - p_recon),[],'all'));
colorbar; 
title       ("Difference Image [dB re maxval]")
subtitle    (100*geomean(abs(p_recon2 - p_recon)./abs(p_recon),'all') + "% Error (geometric mean)")

nexttile; 
imagesc     (abs(p_recon2 - p_recon));
title       ("Difference Image [timeseries units]")

            %% Back-Calculate the required pressure-field

M2          = zeros(size(rf_double));
M2          (cs_inds.time,cs_inds.space) = RF2;

p_backcalc  = real(fftshift(ifftn(ifftshift(M2))));
p_backcalc  = p_backcalc(data.frame_size+1:end,:);

figure      ("Position",[800 10 800 1300])
tiledlayout (3,1); 
nexttile;
imagesc     (data.rfdata); 
title       ("Original Pressure Field");

nexttile;
imagesc     (p_backcalc); 
title       ("Back Calculated Pressure Field");

nexttile;
imagesc     ((abs(p_backcalc - data.rfdata))); 
colorbar; colormap(flip(pucolors.purplebone))
title       ("Difference Image")

            %% Helper Functions
function fig_out = db_plot(k_spc_img)

    img_in_db = 20.*log10(abs(k_spc_img));
    fig_out = imagesc(img_in_db);
    fig_out.Parent.YDir = 'normal';
    daspect([1 1 1]);
    colormap(pucolors.cmrmap)
    colorbar

end

function plot_two_leaves(leaf1, leaf2, n_zero_pix)
    figure("Position",[600 600 1300 600]); 
        tiledlayout(1,2); nexttile;
        imagesc(leaf1); colormap(pucolors.purplebone); daspect([1 1 1]);
        nexttile;
        imagesc(leaf2); colormap(pucolors.purplebone); daspect([1 1 1]);
        title("Compression Ratio: " + n_zero_pix/numel(leaf2)*100 + "%")
        subtitle("Zero Pixels: " + n_zero_pix)
end

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