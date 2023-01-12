%Sensor Configuration

classdef ImageArea

properties
    x_arr
    y_arr
    pulse_ind
end

methods

    %Constructor
    function obj = ImageArea(xmax, xmin, xint, ymax, ymin, yint)
    
        obj.x_arr = (xmin:xint:xmax).';
        obj.y_arr = (ymin:yint:ymax).';

    end

    %GPU version of imaging area object.
    function obj = gpuImageArea(obj)

        obj.x_arr = gpuArray(obj.x_arr);
        obj.y_arr = gpuArray(obj.y_arr);
        obj.pulse_ind = gpuArray(obj.pulse_ind);
    
    end
end



end