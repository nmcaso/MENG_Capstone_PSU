classdef SensorArray
%Sensor Array object

properties
    x0      % x-coordinates of the transducers
    z0      % y-coordinates of the transducers (with naming convention of the y axis being depth)
end

methods
    
    function obj = SensorArray(filepath)
    %Constructor - load the Sensor Array config

        % Input Checks
        switch class(filepath)
            case 'string'
                load(filepath,'z0','x0')
            case 'char'
                load(string(filepath),'z0','x0')
            otherwise
            warning("Error: file path is not a string that points to a dataset .mat file")
        end

        % Convert from mm to meters
        obj.x0  = x0*1e-3;
        obj.z0  = z0*1e-3;

        % Correct faulty sensors
        faulty = abs(diff(obj.x0)) > 8*mean(obj.x0);
        if any(faulty)
            %interpolate between the adjacent two sensors if there are
            %faulty sensor positions
            faultyinds = find(faulty);
            obj.x0(faultyinds(2:2:end)) = (obj.x0(faultyinds(1:2:end)) + obj.x0(faultyinds(2:2:end)+1))/2;
            obj.z0(faultyinds(2:2:end)) = (obj.z0(faultyinds(1:2:end)) + obj.z0(faultyinds(2:2:end)+1))/2;
        end
    end

    function obj = gpuSensorArray(obj)
    %Convert to GPU arrays if desired
        obj.x0  = gpuArray(obj.x0);
        obj.z0  = gpuArray(obj.z0);
    end
end
end
