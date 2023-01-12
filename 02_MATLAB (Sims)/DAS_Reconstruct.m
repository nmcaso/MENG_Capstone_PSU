            %% Setup
            clear; clc
            refdir = cd; refdir = refdir(1:55);
            dirs = refdir + "\" + strtrim(string(ls([refdir '\*0*'])));
            for ii = dirs'
            addpath(ii);
            end

            workdir = cd;
            dirs = strtrim(string(ls([workdir '\*0*'])));
            for ii = 1:length(dirs)
            addpath(string(workdir) + "\" + dirs(ii))
            end            
            %% Dataset and Imaging Area

dta         = Dataset("ModelRun.mat",0); % dta = dta.rfcaster('double');
sensarr     = SensorArray("ModelRun.mat");

if size(dta.rfdata,1) ~= 2125
    dta.rfdata = [dta.rfdata; zeros(2125-size(dta.rfdata,1),512)];
end

mins        = -0.2; intval = 1e-3; maxs = -mins - intval;
area1       = ImageArea(maxs, mins, intval, maxs, mins, intval);
indmat      = IndexMatrix(sensarr,area1,dta,"index",false); disp(indmat.Times.total);

indmat.M(indmat.M <= 0) = 1;
indmat.M(indmat.M > numel(dta.rfdata)) = numel(dta.rfdata);

            %% Run DAS function
[DASIMG, performance]     = DAS_index(dta,indmat);
disp(performance)

DASIMG                    = reshape(DASIMG,length(area1.x_arr),length(area1.y_arr));
                       
           %% Figures

imagesc((DASIMG)); colormap("pucolors.inferno")
title("Function: DAS Index")
colorbar
daspect([1 1 1])
