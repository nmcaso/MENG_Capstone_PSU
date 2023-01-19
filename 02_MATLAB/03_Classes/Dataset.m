% Dataset Class

classdef Dataset

%dataset properties
properties
    fs
    dt
    c
    rfdata
    pulse_ind
end

methods
    
    %Constructor - load the dataset
    function obj = Dataset(filepath,pulse_ind)
        switch class(filepath) 
            case "string"
            load(filepath,'fs','sos','rfdata')
            otherwise
            warning("Error: file path is not a string that points to a dataset .mat file")
        end

        obj.fs  = fs;
        obj.dt  = 1/(fs*1e6);
        obj.c   = sos*1000;
        obj.rfdata = rfdata./size(rfdata,2);
        obj.pulse_ind = pulse_ind;

    end

    %Convert to GPU arrays if desired
    function obj = gpuDataSet(obj)

        obj.fs  = gpuArray(obj.fs);
        obj.dt  = gpuArray(obj.dt);
        obj.c   = gpuArray(obj.c);
        obj.rfdata = gpuArray(obj.rfdata);

    end

    %Change the rfdata data type (in case you need to for sparse mmult)
    function obj = rfcaster(obj,newclass)
        obj.rfdata = cast(obj.rfdata,newclass);
    end

end

end