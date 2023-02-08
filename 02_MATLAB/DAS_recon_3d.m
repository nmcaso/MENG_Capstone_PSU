            %% Setup
            clear; clc; close all;
            for folder = strtrim(string(ls('*0*'))).'
                addpath(folder);
            end

nvidiaDev   = gpuDevice;
data        = Dataset("leaf", 0); 
data        = data.rfcaster('double');
sensarr     = SensorArray("leaf");
interpolt   = false;

mdl.xmin    = -0.02;                    mdl.xres    = 0.0001;
mdl.xmax    = -mdl.xmin - mdl.xres;     mdl.ymin    = mdl.xmin;
mdl.yres    = mdl.xres;                 mdl.ymax    = mdl.xmax;
squarea     = ImageArea(mdl);

ind_mat     = IndexMatrix(sensarr, squarea, data, "index", interpolt);
ind_mat.M   = ind_mat.M - 1;

            %% Function

numavgs     = 200;
img         = zeros(length(squarea.x_arr),length(squarea.y_arr), numavgs);
            for ii = -numavgs/2+1:numavgs/2
img(:,:,ii+numavgs/2)    = dasindex_original(data, squarea, sensarr, ii/10000);
            disp("Iteration " + (numavgs - (ii + numavgs/2)));
            end

            %% Processing
IMG         = fftshift(fftn(fftshift(img),[800 800 400]));
C_IMG       = IMG.*conj(IMG);

corr_img    = abs(ifftshift(ifftn(ifftshift(C_IMG),[400 400 200])));
corr_img    = corr_img./max(corr_img,[],'all');

% corr_img(corr_img < 1e-2) = 0;

Final_img   = abs(corr_img).*img;

            %% Display

view        = viewer3d(BackgroundColor="white", BackgroundGradient=false, Lighting="off");
img3d       = volshow(Final_img, RenderingStyle="GradientOpacity",Parent=view);

