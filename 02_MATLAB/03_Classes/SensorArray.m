classdef SensorArray
%Sensor Array object

properties
    x0      =[] % x-coordinates of the transducers
    z0      =[] % y-coordinates of the transducers (with naming convention of the y axis being depth)
end

properties (Dependent)
    nSensors
    Array
end

methods
    
    function obj = SensorArray(filepath)
    %Constructor - load the Sensor Array config

        % Input Checks
        switch class(filepath)
            case {'string' 'char'}
                load(string(filepath),'x0','z0');
            otherwise
            error("File path is not a string that points to a dataset .mat file");
        end

        % Convert from mm to meters
        obj.x0  = x0(:).'*1e-3;
        if exist("z0",'var')
            obj.z0  = z0(:).'*1e-3;
        end

        % Correct faulty sensors
        faulty = abs(diff(obj.x0)) > 8*mean(obj.x0 + 1) - 1;
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

%get method for dependent properties
methods
    
    function nsens = get.nSensors(obj)
        if min(size(obj.x0)) > 1
            nsens = size(obj.x0);
        else
            nsens = length(obj.x0);
        end
    end

    function array = get.Array(obj)
        
        array(size(obj.x0)) = struct;
        for ii = 1:size(obj.x0,1)
            for jj = 1:size(obj.x0,2)
                array(ii, jj).x = obj.x0(ii, jj);
                array(ii, jj).z = obj.z0(ii, jj);
                array(ii, jj).index = [ii jj];
            end
        end
    
    end

end

end
