            %% A script specifically for making speed comparisons
            clear; clc; close;
            for ii = strtrim(string(ls))'; if isfolder(ii); addpath(ii); end; end

            %% List of the functions

min_speed   = 0.0033;
grow_factor = 1.73;
function_list = {   "DAS_GPU_index" "CUDA_DAS_index" "DAS_original" "dasindex_original" "DAS_index" "dasindex" "DAS_mmult" "DAS_mmult_gpu" "DAS_elementalLoop";
                    "index"         "index"          "index"        "index"             "index"     "index"     "matrix"   "matrix" "index";
                    true true false false false false false true false};

            %% Function Inputs SET TO 256 PIXELS
nvidiaDev   = gpuDevice;
sensarr     = SensorArray("leaf"); sensarr.x0 = double(sensarr.x0(2:2:end)); sensarr.z0 = double(sensarr.z0(2:2:end));
interpolt   = false;

            %% CPU Functions comparison
num_size_inc = 15;
num_avgs    = 35;

results(length(function_list)) = struct; 

function_index = 0;
for func = function_list
    %clean up from the last function
    clear mdl datas
    
    %increment the index
    function_index = function_index + 1;

    %set the model initial image array
    mdl         = struct(               ...
                "xmin", -.5e-2,          ...
                "xres", 1e-4,           ...
                "xmax", .5e-2 - 1e-4,    ...
                "ymin", -.5e-2,          ...
                "yres", 1e-4,           ...
                "ymax", .5e-2 - 1e-4     );
    %load the dataset
    data        = Dataset("leaf",0); data = data.rfcaster('double');
    data.rfdata = data.rfdata(:,2:2:end);
    
    %preallocate some variables
    results(function_index).time_fps = zeros(num_size_inc,1);
    results(function_index).time_pps = zeros(num_size_inc,1);
    results(function_index).img_size = zeros(num_size_inc,1);
    results(function_index).label = func{1};

    slow_flag = false;
    %increment between image sizes
    for size_inc = 1:num_size_inc
        
        %if the last one was really slow (less than half a frame per second)
        if slow_flag
            continue %try again
        end
        
        %clean up from the last size
        clear imat DASIMG args

        %display progress messagecl
        fprintf("Now Processing " + func{1} + " size increment " + size_inc + " of " + num_size_inc + "...");
        
        %initialize the imagearea and index matrix
        squarea = ImageArea(mdl);
        imat = DelayMatrix(sensarr, squarea, data, func{2}, 0);

        if func{3} %send the stuff that's supposed to be on the GPU to the GPU
args    = argument_parse(func{1}, data.gpuDataSet, imat.Send2GPU, squarea, sensarr, 0, nvidiaDev);
        else
args    = argument_parse(func{1}, data, imat, squarea, sensarr, 0, nvidiaDev);
        end

        %set the input arguments
        tic
        for ii      = 1:num_avgs
            DASIMG = feval(func{1}, args{:});
        end
        time_t = toc;
        fprintf(1./(time_t/num_avgs) + " FPS\n");
    
        %calculate and save results
        results(function_index).time_fps(size_inc) = 1./(time_t/num_avgs);
        results(function_index).img_size(size_inc) = numel(DASIMG);
        results(function_index).time_pps(size_inc) = numel(DASIMG)./time_t;
        disp(sqrt(numel(DASIMG)));

        %increase the model size by 5%
        for mdlfield = fields(mdl)'
            if all(mdlfield{:} ~= 'xres') || all(mdlfield{:} ~= 'yres') 
                mdl.(mdlfield{:}) = mdl.(mdlfield{:}) * (sqrt(grow_factor)/2-1/2 + 1);
            end
        end
    
        slow_flag = 1./(time_t/num_avgs) < min_speed;
       
    end
end

save("05_DataOut\Speed_Comparisons256sens.mat","results");

        %% helper function
function out = argument_parse(func_name, dset, imat, area, sens, param, device)

    switch func_name
        case "DAS_original"
            out = {dset, area, sens};
        case "DAS_elementalLoop"
            antiphase(1,1,:) = uint32(2125*(0:size(dset.rfdata,2)-1));
            nophase = imat.M - antiphase;
            out = {dset, nophase};
        case "dasindex_original"
            out = {dset, area, sens, param};
        case "DAS_index"
            out = {dset, imat};
        case "dasindex"
            out = {imat.M - 1, dset.rfdata};
        case "DAS_GPU_index"
            out = {dset, imat};
        case "CUDA_DAS_index"
            out = {imat.M - 1, dset.rfdata};
        case "DAS_mmult"
            out = {dset, imat, area};
        case "DAS_mmult_gpu"
            out = {dset, imat, area, device};
        otherwise
            error("no");
    end

end