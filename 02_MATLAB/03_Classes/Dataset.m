classdef Dataset
% Dataset Class

properties
    fs          % sample rate
    dt          % 1/fs
    c           % speed of sound
    rfdata      % raw time data
    pulse_ind   % pulse index if desired
end

properties (Dependent)
    frame_size
    n_frames
end

methods
    
    function obj = Dataset(filepath, pulse_ind, type)
    %Constructor - load the dataset

        arguments
            filepath string {mustBeFile}
            pulse_ind (1,1) {mustBeScalarOrEmpty} = 0
            type string {mustBeMember(type, ["time" "frequency"])} = "time"
        end

        load(filepath,'fs','sos','rfdata') %load the file
        
        if ~exist('fs','var'); error("Data file missing sample rate fs"); end
        if ~exist('sos','var'); error("Data file missing speed of sound sos"); end
        if ~exist('rfdata', 'var'); error("Data file is missing rfdata matrix"); end

        
        obj.fs      = fs;
        obj.dt      = 1/(fs*1e6);
        obj.c       = sos*1000;
        obj.rfdata  = rfdata./size(rfdata,2);
        obj.pulse_ind = pulse_ind;
        
    end

    function obj = gpuDataSet(obj)
    %Convert to GPU arrays if desired
        obj.fs  = gpuArray(obj.fs);
        obj.dt  = gpuArray(obj.dt);
        obj.c   = gpuArray(obj.c);
        obj.rfdata = gpuArray(obj.rfdata);
    end

    function obj = rfcaster(obj,newclass)
    %Change the rfdata data type (in case you need to for sparse mmult)
        obj.rfdata = cast(obj.rfdata,newclass);
    end

end

%set/get methods for properties
methods
    
    function obj = set.fs(obj, inp)

    function frsize = get.frame_size(obj)
        frsize = size(obj.rfdata,1);
    end

    function nfr = get.n_frames(obj)
        if numel(size(obj.rfdata)) == 2
            nfr = 1;
        elseif numel(size(obj.rfdata)) == 3
            nfr = size(obj.rfdata,3);
        end
    end

end

end