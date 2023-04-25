            %% Interpolated DAS
            clear;clc;close all;
            for fldr = strtrim(string(ls))'; if isfolder(fldr); addpath(fldr); end; end

            %% Setup
data        = Dataset("crossinghair.mat"); data = data.rfcaster("double");
sens        = SensorArray("crossinghair.mat"); 
sens.x0     = cast(sens.x0,"double");
sens.z0     = cast(sens.z0,"double");
img         = ImageArea(0.0199, -0.02, 1e-4, 0.0199, -0.02, 1e-4);
img_recon0  = dasindex_original(data, img, sens, 0);

factor      = 2;
RF          = fft2(data.rfdata, size(data.rfdata,1)*factor, size(data.rfdata,2)*factor);
data.rfdata = real(ifft2(RF));

sens.z0     = real(ifft(fft(sens.z0, factor*length(sens.z0))));
sens.x0     = real(ifft(fft(sens.x0, factor*length(sens.x0))));

            %% Run DAS!
img_recon   = dasindex_original(data, img, sens, 0);

tiledlayout(1,3);
            nexttile; imagesc(img_recon0); colorbar; title("Original")
            nexttile; imagesc(img_recon); colorbar; title("Interp")
            nexttile; imagesc(abs(img_recon) - abs(img_recon0)); colorbar; title("diff")