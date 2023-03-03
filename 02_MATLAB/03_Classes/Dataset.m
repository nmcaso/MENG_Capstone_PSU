classdef Dataset
% Dataset Class

properties
    fs          % sample rate
    dt          % 1/fs
    c           % speed of sound
    rfdata      % raw time data
    pulse_ind   % pulse index if desired
end

methods
    
    function obj = Dataset(filepath,pulse_ind)
    %Constructor - load the dataset

        switch class(filepath) 
            case "string"
            load(filepath,'fs','sos','rfdata') %load the file
            case 'char'
            load(string(filepath),'fs','sos','rfdata')
            otherwise
            warning("Error: file path must lead to a valid dataset .mat file")
        end
        
        try %import the variables from the data file
        obj.fs      = fs;
        obj.dt      = 1/(fs*1e6);
        obj.c       = sos*1000;
        obj.rfdata  = rfdata./size(rfdata,2);
        obj.pulse_ind = pulse_ind;
        catch
        error("Something was wrong with one of the variables. Please check the dataset file.")
        end
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

end