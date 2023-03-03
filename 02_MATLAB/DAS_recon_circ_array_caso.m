            %% Setup
            clear; clc; close all;
            for ii = strtrim(string(ls))'; if isfolder(ii); addpath(ii); end; end

            % Select the things to run:
c_img           = false;
mex_original    = false;
framerotate     = false;
framerotate_gpu = false;
ind_mat_vec     = false;
ind_mat_element = false;
ind_mat_halfvec = false;
ind_c           = false;
ind_gpu_mat     = false;
ind_gpu_cuda    = false;
mmt_matlab      = false;
mmt_gpu_matlab  = false;

compute_averages= false;
nicerplots      = false;
frotate_comparison = true;

compile_mex_or  = false;
compile_c_mex   = false;
compile_cuda_mex= true;

if nicerplots; ind_c = true; end
if frotate_comparison, ind_c = true; framerotate_gpu = true; end

switches        = [ind_mat_vec ind_mat_halfvec ind_mat_element ind_c ind_gpu_mat ind_gpu_cuda ...
                   mmt_matlab mmt_gpu_matlab c_img mex_original framerotate framerotate_gpu];
gpu_switches    = [ind_gpu_mat ind_gpu_cuda mmt_gpu_matlab framerotate_gpu];
no_ind_matrix   = [c_img mex_original framerotate framerotate_gpu];

numimgs         = sum(switches);

            %% Compile MEX function(s) if desired
            if compile_c_mex
folder      = pwd;
copyfile([folder(1:46) '03_C++\MEX Functions\DAS_Index.c'])
             mex "DAS_Index.c" -I"C:\Program Files\MATLAB\R2022b\extern\include" COPTIMFLAGS="-fwrapv" -output dasindex
            end
            if compile_mex_or
folder      = pwd;
copyfile([folder(1:46) '03_C++\MEX Functions\DAS_Index_Original.cpp'])
            mex "DAS_Index_Original.cpp" -I"C:\Program Files\MATLAB\R2022b\extern\include" COPTIMFLAGS="-O3 -fwrapv" -output dasindex_original
            end
            if compile_cuda_mex
folder      = pwd;
copyfile([folder(1:46) '03_C++\MEX Functions\gpudasindex.cu'])
            mexcuda "gpudasindex.cu" COPTIMFLAGS="-O3 -fwrapv" -output CUDA_DAS_index
            end

            %% Create Dataset, imaging area, and delay matrix
nvidiaDev   = gpuDevice;
data        = Dataset("leaf",0); 
data        = data.rfcaster('double');
sensarr     = SensorArray("leaf"); sensarr.x0 = double(sensarr.x0); sensarr.z0 = double(sensarr.z0);
interpolt   = false;

mdl.xmin    = -0.02;                    mdl.xres    = 0.00006;
mdl.xmax    = -mdl.xmin - mdl.xres;     mdl.ymin    = mdl.xmin;
mdl.yres    = mdl.xres;                 mdl.ymax    = mdl.xmax;
% mdl.zmin    = -0.01;                    mdl.zres    = 0.001;
% mdl.zmax    = -mdl.zmin - mdl.zres;

squarea     = ImageArea(mdl);

            if any(switches(1:8))
ind_mat     = DelayMatrix(sensarr, squarea, data, "index", interpolt); 
            disp("Time to make Index Matrix: " + ind_mat.Times.total)
            end

            if any(gpu_switches) || compute_averages
gpudata     = data.gpuDataSet;
gind_mat    = ind_mat.Send2GPU;
            end
            if mmt_matlab || compute_averages
mult_mat    = DelayMatrix(sensarr,squarea,data,"matrix", interpolt); disp("Time to make MMT Matrix: " + mult_mat.Times.total);
            end
            if mmt_gpu_matlab || compute_averages
gmul_mat    = mult_mat.Send2GPU;
            end

            %% Run the original DAS function without optimization
            if c_img
            tic
MAT_ORIGIN  = DAS_original(data, squarea, sensarr);
            time_orig = toc;
            disp("Original function time: " + time_orig)
colormax    = max(abs(MAT_ORIGIN),[],'all');

            nexttile; imagesc(MAT_ORIGIN); colormap(pucolors.seismic); clim([-colormax colormax]) 
            end

            %% DAS original function compiled
            if mex_original
            tic
MAT_MEX_ORI = dasindex_original(data, squarea, sensarr, 0);
            time_compiled_orig = toc;
            disp("Compiled Original DAS function: " + time_compiled_orig);
colormax    = max(abs(MAT_MEX_ORI),[],'all');

            nexttile; imagesc(MAT_MEX_ORI); colormap(pucolors.seismic);
            end

            %% Frame rotate on CPU
            if framerotate
            tic
MAT_ROTATE  = DAS_frame_rotate(data, squarea, sensarr, 9.1, true);
            time_frotate = toc;
            disp("Frame Rotation time: " + time_frotate)

            nexttile; imagesc(abs(MAT_ROTATE)); colormap(pucolors.purplebone); 
            end
            %% Frame rotate on GPU
            if framerotate_gpu
            gpudata.rfdata = gpuArray(single(gpudata.rfdata));
            tic
MAT_GPU_ROT = DAS_frame_rotate(gpudata, squarea, sensarr, 9.1, false);
            time_gpu_frotate = toc;
            disp("GPU Frame Rotation time: " + time_gpu_frotate)
            colormax = gather(max(abs(MAT_GPU_ROT),[],'all'));

            nexttile; imagesc(MAT_GPU_ROT); colormap(pucolors.seismic); clim([-colormax colormax]) 
            end
            %% Run indexing function normally on CPU
            if ind_mat_vec
            tic;
MAT_VEC_IMG = DAS_index(data,ind_mat);
            time_matlab_vec = toc;
            disp("DAS Indexing (Matlab Fully Vectorized) time: " + time_matlab_vec)

            switch ind_mat.type
                case "Index2D"
                
                nexttile; imagesc(abs(MAT_VEC_IMG)); colormap(pucolors.purplebone);
            
                case "Index3D"
                
                max_intens = max(abs(MAT_VEC_IMG),[],'all'); 
                num_intspts = 10;
                intens_scale = max_intens/num_intspts:max_intens/num_intspts:max_intens;
                alpha = (0:num_intspts-1)./num_intspts;
                alpha = 100.^(alpha);
                alpha = alpha./max(alpha);

                queryPoints = linspace(min(intens_scale),max(intens_scale),256);
                alphamap = interp1(intens_scale,alpha,queryPoints)';


view        = viewer3d(BackgroundColor="white", BackgroundGradient=false, Lighting="off");
            volshow(abs(MAT_VEC_IMG), Alphamap=alphamap, Parent=view);

            end
            
            end
            
            %% Run half-vectorized function on CPU
            if ind_mat_halfvec
antiphase(1,1,:) = uint32(2125*(0:511));
non_phased_M = ind_mat.M - antiphase;
            tic;
MAT_HALF_IMG= DAS_half_vec_loop(data, non_phased_M);
            time_matlab_halfvec = toc;
            disp("DAS Indexing (Matlab Partially Vectorized) time: " + time_matlab_halfvec);
            
            nexttile; imagesc(abs(MAT_HALF_IMG)); colormap(pucolors.purplebone); 
            end

            %% Run completely elemental matlab function on CPU
            if ind_mat_element
antiphase(1,1,:) = uint32(2125*(0:511));
non_phased_M = ind_mat.M - antiphase;
            tic;
MAT_ELE_IMG = DAS_elementalLoop(data, non_phased_M);
            time_matlab_elem = toc;
            disp("DAS Indexing (Matlab Elemental) time: " + time_matlab_elem);
            
            nexttile; imagesc(abs(MAT_ELE_IMG)); colormap(pucolors.purplebone); 
            end

            %% Run indexing via C MEX function:
            if ind_c
mat         = ind_mat.M - 1; %because C is zero-indexed
            
            tic
C_IMAGE     = dasindex(mat, data.rfdata);
            time_c = toc; disp("DAS Index in C: " + time_c)

nexttile; imagesc(abs(C_IMAGE)); colormap(pucolors.cvidis);
            end

            %% Run indexing function on GPU
            if ind_gpu_mat
MAT_GPU_IMG = DAS_GPU_index(gpudata, gind_mat);
            nvidiaDev.wait();
            time_gpu = toc;
            disp("GPU (Matlab) performance: " + time_gpu)

            nexttile; imagesc(abs(MAT_GPU_IMG)); colormap(pucolors.cvidis);
            end

            %% Run indexing function on CUDA GPU
            if ind_gpu_cuda
smat        = gind_mat.M-1; % again, because C is zero-indexed
            nvidiaDev.wait(); %wait for the subtraction to end on all threads

            tic
CUDA_IMG    = CUDA_DAS_index(smat, gpudata.rfdata);
            time_cuda = toc; disp("GPU (CUDA) performance: " + time_cuda);
            % no nvidiaDev.wait() function necessary here; CUDA_DAS_index
            % includes the cuda device synchronization calrfdata      = mxGPUCreateFromMxArray(prhs[1]);ls.

            nexttile; imagesc(abs(CUDA_IMG)); colormap(pucolors.cvidis);
            end

            %% Run mmult DAS function normally on CPU
            if mmt_matlab
            tic;
MAT_MUL_IMG = DAS_mmult(data, mult_mat, squarea);
            time_mmult_cpu = toc;
            disp("Matrix Mult DAS Time: " + time_mmult_cpu);   

            nexttile; imagesc(abs(MAT_MUL_IMG)); colormap(pucolors.cvidis);
            end

            %% Run mmult DAS function on GPU
            if mmt_gpu_matlab
            tic
MAT_MUL_GPU = DAS_mmult(gpudata,gmul_mat,squarea);
            nvidiaDev.wait();
            time_mmult_gpu = toc;
            disp("GPU MMULT Time: " + time_mmult_gpu)

            nexttile; imagesc(abs(MAT_MUL_GPU)); colormap(pucolors.cvidis);
            end

            %% Averages
            if compute_averages
numavgs     = 100;
cpumnd = zeros(numavgs,1); halfvec = cpumnd; cudaind = cpumnd;
gpuind = cpumnd; cpumult = cpumnd; gpumult = cpumnd; cpucnd = cpumnd; elemt = cpumnd;
            
            for aa = 1:numavgs %Normal indexing function on CPU
            tic;
IMG1        = DAS_index(data, ind_mat);
cpumnd(aa)  = toc;
            end; disp("Indexing CPU Completed")

antiphase(1,1,:) = uint32(2125*(0:511));
non_phased_M = ind_mat.M - antiphase;
%             for aa = 1:numavgs %Elemental Loop Indexing
%             tic;
% IMG2        = DAS_elementalLoop(data, non_phased_M);
% elemt(aa)   = toc;
%             end; disp("Elemental Indexing Complete")

            for aa = 1:numavgs %Partially vectorized loop indexing
            tic;
IMG3        = DAS_half_vec_loop(data, non_phased_M);
halfvec(aa) = toc;
            end; disp("Half-Vectorized indexing Complete")

mat         = ind_mat.M - 1;
            for aa = 1:numavgs %Normal indexing function with MEX
            tic
IMG4        = dasindex(mat, data.rfdata);
            ccputime = toc;
cpucnd(aa)  = ccputime;
            end; disp("C MEX Indexing CPU Completed")

            for aa = 1:numavgs %Normal indexing function on GPU
            tic
IMG5        = DAS_GPU_index(gpudata, gind_mat);
            nvidiaDev.wait();
gpuind(aa)  = toc;
            end; disp("Indexing GPU Completed")

smat        = gind_mat.M - 1; nvidiaDev.wait();
            for aa = 1:numavgs
            tic
IMG6        = CUDA_DAS_index(smat, gpudata.rfdata);
cudaind(aa) = toc;
            end; disp("CUDA Indexing Complete");

            for aa = 1:numavgs %Sparse Matrix Multiplication on CPU
            tic;
IMG7        = DAS_mmult(data, mult_mat, squarea);
cpumult(aa) = toc;
            end; disp("MMULT CPU Completed")

            for aa = 1:numavgs %Sparse Matrix Multiplication on GPU
            tic;
IMG8        = DAS_mmult(gpudata,gmul_mat,squarea);
            nvidiaDev.wait();
gpumult(aa) = toc;
            end; disp("MMULT GPU Completed\n")

trunc       = 1;
avg_iele    = sum(elemt(trunc:end))/(numavgs-trunc+1);
avg_icpu    = sum(cpumnd(trunc:end))/(numavgs-trunc+1);
avg_half    = sum(halfvec(trunc:end))/(numavgs-trunc+1);
avg_icpuc   = sum(cpucnd(trunc:end))/(numavgs-trunc+1);
avg_igpu    = sum(gpuind(trunc:end))/(numavgs-trunc+1);
avg_cuda    = sum(cudaind(trunc:end))/(numavgs-trunc+1);
avg_mcpu    = sum(cpumult(trunc:end))/(numavgs-trunc+1);
avg_mgpu    = sum(gpumult(trunc:end))/(numavgs-trunc+1);

fprintf("Average processing time indexing elemental     %.4f seconds.\n", avg_iele);
fprintf("Average processing time indexing half-vec:     %.4f seconds.\n", avg_half);
fprintf("Average processing time vectorized on CPU:     %.4f seconds.\n", avg_icpu);
fprintf("Average processing time indexing with C MEX:   %.4f seconds.\n", avg_icpuc);
fprintf("Average processing time indexing on GPU:       %.4f seconds.\n", avg_igpu);
fprintf("Average processing time matrix mult CPU:       %.4f seconds.\n", avg_mcpu);
fprintf("Average processing time matrix mult GPU:       %.6f seconds.\n", avg_mgpu);
fprintf("Average processing time CUDA Indexing:         %.6f seconds.\n", avg_cuda);
            end
            %% Figures
            if nicerplots  
            close all
mx1 = max(abs(C_IMAGE(:)));

c_img = abs(C_IMAGE/mx1); %normalize

a = 1;
nonintlog = (20*log10(1+a*c_img)./log10(1+a));

fig1    = image_plot(nonintlog, squarea);
fig2    = dataset_plot(data);

saveas(fig1,"05_DataOut\SampleIMG_" + size(c_img,1) + "x" + size(c_img,2)+".png");
saveas(fig2,"05_DataOut\Sample_rfplot.png")
            end

            %% Frame Rotation Comparison Figure
            if frotate_comparison
            close all

mx1         = max(abs(C_IMAGE(:)));
c_img       = abs(C_IMAGE./mx1); %normalize

mx2         = max(abs(MAT_GPU_ROT(:)));
r_img       = abs(MAT_GPU_ROT'./mx2); %normalize

b = 1;
cimg        = (20*log10(1+b*c_img)./log10(1+b));
rotimg      = (20*log10(1+b*r_img)./log10(1+b));

img1 = image_plot(rotimg, squarea);
img2 = image_plot(cimg, squarea);
img1.Children(2).CLim = [1 20];
img1.Children(2).Subtitle = [];
img2.Children(2).CLim = [1 20];
img2.Children(2).Subtitle = [];

saveas(img1, "05_DataOut\FR_Compare_rotated.png")
saveas(img2, "05_DataOut\FR_Compare_indexed.png")
            end

            %% Functions

function plot_out = dataset_plot(data)

    plot_out    = figure;
    data_dB     = 20*log10(abs(data.rfdata)/1e-6);

    plt         = imagesc(data_dB); colormap(pucolors.viridis);
    clim        ([50 120]);
    ax          = plt.Parent;
    ax          .YDir = 'normal';

    c           = colorbar;
    c.          Label.String = "Sound Pressure Level (dB re 1 \muPa)";    
    
    ax.         YTickLabel = 1e6*(ax.YTick * data.dt);

    xlabel      ("Sensor Number");
    ylabel      ("Time of Arrival, \mu"+"s (1\times 10^-^6 seconds)");
    daspect     ([1 3.5 1]);

end

function plot_out = image_plot(img, area, varargin)
    
    plot_out    = figure('Position', [100 100 500 500]);
    plt         = imagesc(img); colormap(pucolors.cmrmap);
%     clim        ([0.2 2.4]);
    ax          = plt.Parent;
    
    [xsz, ysz]  = size(img);
    
    xtvec       = 0:xsz/10:xsz;
    xtvec       = [1 xtvec(2:end)];
    ax.XTick    = xtvec;
    ax.XTickLabel = area.x_arr(round(xtvec))*1E2;
    
    ytvec       = 0:ysz/10:ysz;
    ytvec       = [1 ytvec(2:end)];
    ax.YTick    = ytvec;
    ax.YTickLabel = area.y_arr(round(ytvec))*1E2;
    
    c           = colorbar;
    c.          Label.String = "Normalized Image Magnitude (dB re 1 Pa)";
    
    xlabel      ("X-position (mm)");
    ylabel      ("Y-position (mm)");
    subtitle    (size(img,1) + " by " + size(img,2) + " Pixels (" + numel(img)/1e6 + " MP)");
    daspect     ([1 1 1]);



end