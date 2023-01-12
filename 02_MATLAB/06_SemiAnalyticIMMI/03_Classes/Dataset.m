% Dataset Class
classdef Dataset

%dataset properties
properties
    fs
    dt
    c
    rfdata
    pulse_ind
    sensLocs
    sensCount
end

methods
    
    %Constructor - load the dataset
    function obj = Dataset(filepath,pulse_ind)
        switch class(filepath) 
            case "string"
            load(filepath,'fs','sos','rfdata','z0','x0')
            otherwise
            warning("Error: file path is not a string that points to a dataset .mat file")
        end

        obj.fs  = fs;
        obj.dt  = 1/(fs*1e6);
        obj.c   = sos*1000;
        obj.rfdata = rfdata;
        obj.pulse_ind = pulse_ind;
        obj.sensLocs = [x0;z0]*1e-3;
        obj.sensCount = size(z0);
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
    
    %Reshape the rfdata to be a 1D vector
    function flatRF = rfFlatten(obj,frame)
        flatRF = reshape(obj.rfdata(:,:,frame),size(obj.rfdata,1)*size(obj.rfdata,2),1);
    end

end

end
