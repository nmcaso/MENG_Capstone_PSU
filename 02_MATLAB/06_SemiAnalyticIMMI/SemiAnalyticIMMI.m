%%
% This program will mimic the methods in the Fast Semi-Analytical paper
% DOI: 10.1109/TMI.2010.2044584
%
% Written 2022 by Paul Klippel as a part of the SIMBA lab at Penn State
        %% Filepaths for the program
    addpath("01_DataIn\"); addpath("02_Functions\"); addpath("03_Classes\");
            addpath("04_DataOut\");
            
        %% Dataset and Imaging region generation
file = "phantom_1mm.mat";
xlim = [-1e-3,1e-3];
ylim = xlim;
res = 1e-4;
global GAM
GAM = 0.1;

imagingData = Dataset("leaf",0); imagingData = imagingData.rfcaster('double');
imageRegion = Image(xlim(1),xlim(2),ylim(1),ylim(2),res);
        %% Interpolate and calculate the forward matrix
g = zeros(size(imagingData.rfdata,1)*size(imagingData.rfdata,2),imageRegion.totalPix);

for sensor = 1:imagingData.sensCount
    x0 = imagingData.sensLocs(1,sensor);
    y0 = imagingData.sensLocs(2,sensor);
    for t_ind = 1:size(imagingData.rfdata,1)
        R = t_ind/imagingData.fs*imagingData.c;
        zMaxLast = 0;
        tSens = t_ind+(size(imagingData.rfdata,1)*(sensor-1));
        for theta = 0:1e-6:2*pi %for now interpolating the whole circle
            NN = Neighbors(R,theta,imageRegion,imagingData.sensLocs(:,sensor));
            if length(NN.xNeighbors)<2 || length(NN.yNeighbors)<2
                %Outside Imaging Region
                continue
            end
            if NN.nnIndex(NN.farthest)==zMaxLast
                %Same triangle
                zMaxLast = NN.nnIndex(NN.farthest);
                continue
            else
                zMaxLast = NN.nnIndex(NN.farthest);
            end
            'here'
            g = TriangleInterp(NN,g,tSens);
            clear NN;
        end
    end
end
        %% Invert the matrix using psuedo-inversion
Mpsuedo = (g'*g)^(-1)*g';

%One frame for image calculation, reshaped to be a 1D vector
pFlat = imagingData.rfFlatten(1);
sol = Mpsuedo*pFlat;

sol2 = g\pFlat; %checking with MATLAB's solution

IMG = imageRegion.FinalShape(sol);
IMG2 = imageRegion.FinalShape(sol2);
figure
imagesc(IMG)
figure
imagesc(IMG2)