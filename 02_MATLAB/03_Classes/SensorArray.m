%Sensor Array object

classdef SensorArray
%Sensor Array properties

properties
    x0
    z0
end

methods
    
    %Constructor - load the Sensor Array config
    function obj = SensorArray(filepath)

        if isa(filepath, "string")
            load(filepath,'z0','x0')
        else
            warning("Error: file path is not a string that points to a dataset .mat file")
        end

        obj.x0  = x0*1e-3;
        obj.z0  = z0*1e-3;

        %detect and fix faulty sensors!
        faulty = abs(diff(obj.x0)) > 8*mean(obj.x0);
        if any(faulty)
            faultyinds = find(faulty);
            obj.x0(faultyinds(2:2:end)) = (obj.x0(faultyinds(1:2:end)) + obj.x0(faultyinds(2:2:end)+1))/2;
            obj.z0(faultyinds(2:2:end)) = (obj.z0(faultyinds(1:2:end)) + obj.z0(faultyinds(2:2:end)+1))/2;
        end

    end

    %Convert to GPU arrays if desired
    function obj = gpuSensorArray(obj)

        obj.x0  = gpuArray(obj.x0);
        obj.z0  = gpuArray(obj.z0);
        
    end

end

end
