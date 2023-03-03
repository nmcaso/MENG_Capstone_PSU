            %% A script to make my figures!
            clear; clc; close all;
            for ii = strtrim(string(ls))'; if isfolder(ii); addpath(ii); end; end
            load("Speed_Comparisons_ele_loop.mat");
results_elem = results(9);
            for ii = fields(results_elem)'
            if length(results_elem.(ii{1})) > 1
results_elem.(ii{1}) = results_elem.(ii{1})(1:12);
            end
            end
results_elem2 = load("Rerun_Bads.mat");

            load("Speed_Comparisons.mat");
results     = [results, results_elem, results_elem2.results(9)];

results(2).label = "bad";

            results     = better_titles(results);
            results     = custom_colors(results);
            results     = correct_pps(results);

tablepath = "C:\Users\cason\OneDrive\Documents\PSU\Project\07_Dissertation Paper\Latex Thesis\00_NCaso_Thesis\Chapter-4\Large_Array_Speed_Comparisons.xlsx";
large_array_table = readtable(tablepath);

mycolors2 = pucolors.magma(5).';

tbl(1).label = "DAS_original";
tbl(1).time_fps = flip(large_array_table.OriginalDAS_non_vectorized_);
tbl(1).img_size = flip(large_array_table.Size_megapixel_*1e6);
tbl(1).time_pps = tbl(1).img_size./tbl(1).time_fps ;
tbl(1).Color = mycolors2(:,2);
tbl(1).Marker = 'x';
tbl(1).Style = '-.';
tbl(1).PlotLabel = "DAS with 2D Delay Matrix";

tbl(2).label = "DAS_fr_cpu";
tbl(2).time_fps = flip(large_array_table.FrameRotationOnCPU);
tbl(2).img_size = flip(large_array_table.Size_megapixel_*1e6);
tbl(2).time_pps = tbl(2).img_size./tbl(2).time_fps;
tbl(2).Color = mycolors2(:,3);
tbl(2).Marker = '>';
tbl(2).Style = ':';
tbl(2).PlotLabel = "Frame Rotation on CPU";

tbl(3).label = "DAS_fr_gpu";
tbl(3).time_fps = flip(large_array_table.FrameRotationOnGPU);
tbl(3).img_size = flip(large_array_table.Size_megapixel_*1e6);
tbl(3).time_pps = tbl(3).img_size./tbl(3).time_fps;
tbl(3).Color = mycolors2(:,4);
tbl(3).Marker = '<';
tbl(3).Style = ':';
tbl(3).PlotLabel = "Frame Rotation on GPU";

loader = load("05_DataOut\Speed_Comparisons256sens.mat");
results256 = loader.results;
results256 = better_titles(results256);
results256 = custom_colors(results256);

clearvars -except results tbl large_array_table results256
cd 05_DataOut\
delete Results*.png
cd ..

            %% Plot the things you want to plot

figs(1)     = plot_fps(results, "DAS_elementalLoop");
figs(2)     = plot_pps(results, "DAS_elementalLoop");
figs(3)     = plot_fps(results, "DAS_original");
figs(4)     = plot_pps(results, "DAS_original");
figs(5)     = plot_fps(results, ["DAS_elementalLoop", "DAS_original","DAS_index"]);
figs(6)     = plot_pps(results, ["DAS_elementalLoop", "DAS_original","DAS_index"]);
figs(7)     = plot_fps(results, ["DAS_original","DAS_index", "DAS_mmult"]);
figs(8)     = plot_pps(results, ["DAS_original","DAS_index", "DAS_mmult"]);
figs(9)     = plot_fps(results, ["DAS_mmult", "DAS_GPU_index", "DAS_mmult_gpu"]);
figs(10)    = plot_pps(results, ["DAS_mmult", "DAS_GPU_index", "DAS_mmult_gpu"]);
figs(11)    = plot_fps(results, ["DAS_original", "dasindex_original", "DAS_index", "dasindex"]);
figs(12)    = plot_pps(results, ["DAS_original", "dasindex_original", "DAS_index", "dasindex"]);
figs(13)    = plot_fps(results, ["DAS_GPU_index", "CUDA_DAS_index", "dasindex"]);
figs(14)    = plot_pps(results, ["DAS_GPU_index", "CUDA_DAS_index", "dasindex"]);
figs(15)    = plot_fps(tbl, ["DAS_original", "DAS_fr_cpu", "DAS_fr_gpu"]);
            figs(15).Children(2).YLabel.String = "Time to Run (seconds)";
figs(16)    = plot_pps(tbl, ["DAS_original", "DAS_fr_cpu", "DAS_fr_gpu"]);
figs(17)    = figure("Name","Speed vs Size",'position', [100,100,400,400]);
            plot(flip(large_array_table.Size_megapixel_),100*flip(large_array_table.PercentageOfMyGPUUsed_),'-r','Marker','o','LineWidth',1.2)
            grid on;
            xlabel("Image Size (megapixels)"); ylabel("Percentage of GPU used")
figs(18)    = figure("Name","Speed Comparison Plot",'position', [100,100,400,400]);
            plot(100*flip(large_array_table.PercentageOfMyGPUUsed_),flip(large_array_table.FrameRotationOnGPU),'-r','Marker','o','LineWidth',1.2)
            grid on;
            ylabel("Imaging Speed (Seconds)"); xlabel("Percentage of GPU used")

figs(19)    = plot_fps(results, [results([1 3:end]).label]);
figs(20)    = plot_pps(results, [results([1 3:end]).label]);
            figs(19).Children(1).Location = 'eastoutside';
            figs(20).Children(1).Location = 'eastoutside';
            figs(20).Children(2).YScale = 'log';
            figs(19).Position = [100 100 800 400];
            figs(20).Position = [100 100 800 400];
figs(21)    = plot_fps(results, [results256.label]);
figs(22)    = plot_pps(results, [results256.label]);
            figs(21).Children(1).Location = 'eastoutside';
            figs(22).Children(1).Location = 'eastoutside';
            figs(22).Children(2).YScale = 'log';
            figs(21).Position = [100 100 800 400];
            figs(22).Position = [100 100 800 400];

            %% Save the results
fno = 1;
for ii = figs
    saveas(ii, "05_DataOut\Results_fig_" + fno,'png');
    fno = fno + 1;
end

            %% Write a table of results
for ii = 1:length(results)
    if ii ~=2
avg_pps(ii) = round(mean(results(ii).time_pps)./1e6,3);
    end
end
avg_pps(2) = avg_pps(10);
avg_pps = avg_pps(1:9);
for ii = 1:length(results256)
avg_pps256(ii) = round(mean(results256(ii).time_pps)./1e6,3);
end

Speed_table = table([results256.PlotLabel]',avg_pps', avg_pps256','VariableNames',{'Method', '512 Transducers Average Speed (MP per second)','256 Transducers Average Speed (MP per second)'});
Speed_table = sortrows(Speed_table,"256 Transducers Average Speed (MP per second)");
            writetable(Speed_table,"05_DataOut\SpeedTable.xlsx");

            %% Helper plot functions
function plot_fig = plot_fps(data, group, varargin)

    if ~isa(data,'struct')
        error("results wasn't a struct!!")
    end
    legend_vector = strings(length(group),1);
    plot_fig = figure("Name","Speed Comparison Plot",'position', [100,100,400,400]);
    hold on;
    color_ind = 1;

    for ii = group(:)'
        plotme = [data.label] == ii;
        lines = plot(data(plotme).img_size./1e6, data(plotme).time_fps, 'Color', data(plotme).Color, 'Marker', data(plotme).Marker, 'LineStyle', data(plotme).Style);
        lines.LineWidth = 1.2;
        legend_vector(color_ind) = data(plotme).PlotLabel;
        color_ind = color_ind + 1;
    end

    ax = plot_fig.Children;
    ax.XScale = 'linear';
    ax.YScale = 'log';
    ax.TickLabelInterpreter = 'latex';
    ax.YTick = 2.^(-2:10);
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.Color = [0.98 0.98 1];
    xlabel("Image Size in Megapixels (1\times 10^6 pixels)")
    ylabel("Frames per Second")

    if color_ind > 2
        legend(legend_vector,'location','NorthOutside')
    end

    fields = varargin(1:2:end);
    vals = varargin(2:2:end);

    for jj = 1:length(fields)
        set(lines,fields{jj},vals{jj})
    end
end

function plot_fig = plot_pps(data, group, varargin)

    if ~isa(data,'struct')
        error("results wasn't a struct!!")
    end
    legend_vector = strings(length(group),1);
    plot_fig = figure("Name","Speed Comparison Plot",'position', [100,100,400,400]);
    hold on;
    color_ind = 1;

    for ii = group
        plotme = [data.label] == ii;
        lines = plot(data(plotme).img_size/1e6, data(plotme).time_pps/1e6, 'Color', data(plotme).Color, 'Marker', data(plotme).Marker, 'LineStyle', data(plotme).Style);
        legend_vector(color_ind) = data(plotme).PlotLabel;
        lines.LineWidth = 1.2;
        color_ind = color_ind + 1;
    end

    ax = plot_fig.Children;
    ax.XScale = 'linear';
    ax.YScale = 'linear';
    ax.TickLabelInterpreter = 'latex';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.Color = [0.98 0.98 1];

    xlabel("Image Size in Megapixels (1\times 10^6 pixels)")
    ylabel("Megapixels per Second")

    if color_ind > 2
        legend(legend_vector,'location','NorthOutside')
    end

    fields = varargin(1:2:end);
    vals = varargin(2:2:end);

    for jj = 1:length(fields)
        set(lines,fields{jj},vals{jj})
    end
end

function in_struct = better_titles(in_struct)

    bad_titles = ["DAS_GPU_index" "CUDA_DAS_index" "DAS_original" "dasindex_original" "DAS_index" "dasindex" "DAS_mmult" "DAS_mmult_gpu" "DAS_elementalLoop" "CUDA_DAS_index"];
    good_titles = ["DAS Delay Matrix on GPU" "DAS Delay Matrix with CUDA" "DAS 2D Delay Matrix" "DAS 2D Delay Matrix compiled in C++" "DAS Index Vectorized" "DAS Index compiled in C" "DAS by Matrix Multiplication" "DAS by Matrix Multiplication on GPU" "DAS Elemental 3D Delay Matrix" "DAS Delay Matrix with CUDA"];

    for ii = 1:length(in_struct)
        for kk = 1:length(bad_titles)
            if in_struct(ii).label == bad_titles(kk)
    in_struct(ii).PlotLabel =  good_titles(kk);
            end
        end
    end
end

function in_struct = custom_colors(in_struct)
    markers = ["o" "+" "x" "square" "diamond" "^" "pentagram" "hexagram" "<", "+"];
    styles = ["--" "--" "-." "-." "-" "-" ":" ":" "-." "--"];
    colors = pucolors.cmrmap(length(in_struct)+4);
    for ii = 1:length(in_struct)
        in_struct(ii).Color = colors(ii+1,:);
        in_struct(ii).Marker = markers(ii);
        in_struct(ii).Style = styles(ii);
    end
end

function in_struct = correct_pps(in_struct)
    for ii = 1:length(in_struct)-1
        in_struct(ii).time_pps = 1./in_struct(ii).time_pps;
    end
end


