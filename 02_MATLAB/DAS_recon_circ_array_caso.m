            %% Setup
            clear; clc;
            for folder = strtrim(string(ls('*0*'))).'
                addpath(folder);
            end

            % Select the things to run:
ind_matlab      = false;
ind_c           = true;
ind_gpu_matlab  = false;
ind_gpu_cuda    = true;
mmt_matlab      = false;
mmt_gpu_matlab  = false;
compute_averages= false;
nicerplots      = false;
cpp_compare     = false;

compile_c_mex   = false;
compile_cuda_mex= true;

numimgs         = ind_matlab + ind_c + ind_gpu_matlab + ind_gpu_cuda ...
                    + mmt_matlab + mmt_gpu_matlab;

%% Compile MEX function(s) if desired

            if compile_c_mex
folder      = pwd;
copyfile([folder(1:46) '03_C++\MEX Functions\DAS_Index.c'])
             mex "DAS_Index.c" -I"C:\Program Files\MATLAB\R2022b\extern\include" COPTIMFLAGS="-fwrapv" -output dasindex
            end
            if compile_cuda_mex
folder      = pwd;
copyfile([folder(1:46) '03_C++\MEX Functions\gpudasindex.cu'])
            mexcuda "gpudasindex.cu" COPTIMFLAGS="-O3 -fwrapv"
            end

figure; tiledlayout(1, numimgs);
            %% Dataset and Imaging Area
nvidiaDev   = gpuDevice;
leafdata    = Dataset("crossinghair",0); leafdata = leafdata.rfcaster('double');
sensarr     = SensorArray("crossinghair");
interpolt   = false;

mins        = -0.02; intval = 0.0001; maxs = -mins-intval;
area1       = ImageArea(maxs, mins, intval, maxs, mins, intval);
ind_mat     = IndexMatrix(sensarr,area1,leafdata, "index", interpolt); disp("Time to make Index Matrix: " + ind_mat.Times.total)

            if ind_gpu_cuda || ind_gpu_matlab || mmt_gpu_matlab || compute_averages
gind_mat    = ind_mat.Send2GPU;
gleafdata   = leafdata.gpuDataSet;
            end
            if mmt_matlab || compute_averages
mult_mat    = IndexMatrix(sensarr,area1,leafdata,"matrix", interpolt); disp("Time to make MMT Matrix: " + mult_mat.Times.total);
            end
            if mmt_gpu_matlab || compute_averages
gmul_mat    = mult_mat.Send2GPU;
            end

            %% Run indexing function normally on CPU
            if ind_matlab
[DASIMG, timemmtcpu]= DAS_index(leafdata,ind_mat);
            disp("DAS Indexing (Matlab) time: " + timemmtcpu)

            nexttile; imagesc(abs(DASIMG)); colormap(pucolors.purplebone); 
            
            end
            %% Run indexing via C function:
            if ind_c
mat         = ind_mat.M - 1;
            
            tic
CDASIMG     = dasindex(mat, leafdata.rfdata);
            timec = toc; disp("DAS Index in C: " + timec)

nexttile; imagesc(abs(CDASIMG)); colormap(pucolors.cvidis);
            end
            %% Run indexing function on GPU
            if ind_gpu_matlab
[G_DASIMG, timegpu]    = DAS_GP(gleafdata,gind_mat);
            disp("GPU (Matlab) performance: " + timegpu)

nexttile; imagesc(abs(G_DASIMG)); colormap(pucolors.cvidis);
            end

            %% Run indexing function on CUDA GPU
smat        = gind_mat.M-1;
nvidiaDev.wait();

            if ind_gpu_cuda

            tic 
CUDA_DASIMG = gpudasindex(smat, gleafdata.rfdata, 2^3, 2^3, 2^4);
            timecuda = toc; disp("GPU (CUDA) performance: " + timecuda);

nexttile; imagesc(abs(CUDA_DASIMG)); colormap(pucolors.cvidis);
            end

            %% Run mmult DAS function normally on CPU
            if mmt_matlab
[M_DASIMG, timemmtcpu]     = DAS_mmult(leafdata,mult_mat,area1); 
            disp("Matrix Mult DAS Time: " + timemmtcpu);   

nexttile; imagesc(abs(M_DASIMG)); colormap(pucolors.cvidis);
            end

            %% Run mmult DAS function on GPU
            if mmt_gpu_matlab
            tic
[GM_DASIMG]    = DAS_mmult(gleafdata,gmul_mat,area1);
            nvidiaDev.wait();
            timemmtcpu = toc;
            disp("GPU MMULT Time: " + timemmtcpu)

nexttile; imagesc(abs(GM_DASIMG)); colormap(pucolors.cvidis);
            end

            %% Averages
            if compute_averages
numavgs     = 50; cpumnd = zeros(numavgs,1); 
gpuind = cpumnd; cpumult = cpumnd; gpumult = cpumnd; cpucnd = cpumnd;
            for aa = 1:numavgs %Normal indexing function on CPU
[DASIMG, timemmtcpu]       = DAS_index(leafdata,ind_mat);
cpumnd(aa)  = timemmtcpu;
            end; fprintf("Indexing CPU Completed\n")
mat         = ind_mat.M - 1;
            for aa = 1:numavgs %Normal indexing function with MEX
            tic
CUDADASIMG  = dasindex(mat, leafdata.rfdata);
            ccputime = toc;
cpucnd(aa)  = ccputime;
            end; fprintf("C MEX Indexing CPU Completed\n")
            for aa = 1:numavgs %Normal indexing function on GPU
[G_DASIMG, timegpu]    = DAS_GP(gleafdata,gind_mat);
gpuind(aa)  = timegpu;
            end; fprintf("Indexing GPU Completed\n")
            for aa = 1:numavgs %Sparse Matrix Multiplication on CPU
[MDASIMG, mperformance]     = DAS_mmult(leafdata,mult_mat,area1);
cpumult(aa) = mperformance;
            end; fprintf("MMULT CPU Completed\n")
            for aa = 1:numavgs %Sparse Matrix Multiplication on GPU
[GMDASIMG, gmperformance]   = DAS_mmult(gleafdata,gmul_mat,area1);
gpumult(aa) = gmperformance;
            end; fprintf("MMULT GPU Completed\n")

avg_icpu    = sum(cpumnd)/numavgs;
avg_icpuc   = sum(cpucnd)/numavgs;
avg_igpu    = sum(gpuind)/numavgs;
avg_mcpu    = sum(cpumult)/numavgs;  
avg_mgpu    = sum(gpumult)/numavgs;

fprintf("Average processing time indexing on CPU: %.4f seconds.\n", avg_icpu);
fprintf("Average processing time indexing with MEX: %.4f seconds.\n", avg_icpuc);
fprintf("Average processing time indexing on GPU: %.4f seconds.\n", avg_igpu);
fprintf("Average processing time matrix mult CPU: %.4f seconds.\n", avg_mcpu);
fprintf("Average processing time matrix mult GPU: %.4f seconds.\n", avg_mgpu);
fprintf("GPU is %.2f times faster than CPU.\n", avg_icpu/avg_igpu);
fprintf("Mmult GPU is %.2f times faster than indexing GPU, and %.2f times faster than indexing with CPU.\n", avg_igpu/avg_mgpu, avg_icpu/avg_mgpu);
            end
            %% Figures
            if nicerplots
mx1 = max(max(DASIMG));
mx2 = max(max(M_DASIMG));

noninterp = abs(DASIMG/mx1); %normalize
interpmm  = abs(M_DASIMG/mx2);

a = 3;
nonintlog = mx1*(log10(1+a*noninterp)./log10(1+a));
yesintlog = mx2*(log10(1+a*interpmm)./log10(1+a));

figure; subplot(1,2,1)
imagesc(nonintlog)
subtitle("Indexed, not Interpolated")
daspect([1 1 1])
subplot(1,2,2)
imagesc(yesintlog)
subtitle("Matrix Mult, Interpolated")
daspect([1 1 1])
            end
            %% Compare to C++ Results
            if cpp_compare
            addpath("C:\Users\cason\OneDrive\Documents\PSU\Project\03_C++\");
a1          = fopen("cpp_indexCpuOut.bin");
a2          = fopen("cpp_indexGpuOut.bin");
a3          = fopen("cpp_spmltCpuOut.bin");
            
c_indexCpu  = fread(a1,160000,'double');
c_indexGpu  = fread(a2,160000*512,'*int32');
c_spmltCpu  = fread(a3,160000,'double');
            
            fclose(a1); fclose(a2); fclose(a3);

dim1        = 400;
dim2        = 400;
dim3        = 1;

C_Index_CPU = reshape(c_indexCpu,[dim1 dim2 dim3]);
C_Mmult_CPU = reshape(c_spmltCpu,[dim1 dim2 dim3]);
try
    C_Index_GPU = reshape(c_indexGpu,[dim1 dim2 dim3]);
catch
    dim3 = 512;
    C_Index_GPU = reshape(c_indexGpu,[dim1 dim2 dim3]);
end


figure('Position',[100 100 1200 500])
tiledlayout(1,3)
nexttile
    imagesc(abs(C_Index_CPU.')); colormap('pucolors.inferno'); title("C++ Index CPU"); colorbar;
    daspect([1 1 1])
nexttile
    try
    imagesc(abs(C_Index_GPU.')); colormap('pucolors.inferno'); title("C++ Index GPU"); colorbar;
    catch
    imagesc(C_Index_GPU(:,:,1)); colormap('pucolors.inferno'); title("C++ Index GPU"); colorbar;
    end
    daspect([1 1 1]);
% nexttile
%     imagesc(abs(DASIMG.')); colormap('pucolors.inferno'); title("MATLAB Index CPU"); colorbar;
%     daspect([1 1 1]);
nexttile
    imagesc(abs(C_Mmult_CPU.')); colormap('pucolors.inferno'); title("C++ Sparse Matrix Multiplication"); colorbar;
    daspect([1 1 1]);
            end

%             %% Backtrack from kspace
% KDASIMG     = fftshift(fftn(fftshift(DASIMG)));
% 
% [xmesh, ymesh] = meshgrid(area1.x_arr, area1.y_arr);
% [thmesh, rmesh] = cart2pol(xmesh, ymesh);
% 
% figure; tiledlayout(1,3); nexttile;
% a = pcolor(thmesh, rmesh, 20*log10(abs(KDASIMG)));
% a.EdgeColor = 'interp';
% 
% leafpad     = [flip(leafdata.rfdata);leafdata.rfdata];
% leafk       = fftn(leafpad);
% 
% nexttile;
% c = imagesc(20*log10(abs(KDASIMG)));
% 
% nexttile;
% b = pcolor(20*log10(abs(leafk)));
% b.EdgeColor = 'interp';
% 
% colormap(pucolors.cvidis)

            %% muck around with indices

% imagesc(CUDA_DASIMG); colormap(pucolors.cvidis)