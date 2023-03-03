classdef ImageArea
%Imaging Region class for a virtual image location

properties
    x_arr       % the x positions
    y_arr       % the y positions
    z_arr       % (optional) depth positions
end

methods

    function obj = ImageArea(xmax, xmin, xres, ymax, ymin, yres, zmax, zmin, zres)
    %Constructor

        switch nargin

            case 1 % the user may input the fields in a single structure
                model = xmax;
                if ~isa(model,'struct')
                    error("The single input option for the class constructor requires a structure.")
                end

                input_fields = string(fields(orderfields(model)));
                predefined_fields = ["xmax"; "xmin"; "xres"; "ymax"; "ymin"; "yres"; "zmax"; "zmin"; "zres"];
                
                if mod(length(input_fields),3) ~= 0 || ~all(input_fields == predefined_fields(1:length(input_fields)))
                    error("Bad input arguments. The structure fields must include: xmin, xres, xmax, ymin, yres, ymax, and optionally zmin, zres, zmax");
                end

                obj.x_arr = double(model.xmin:model.xres:model.xmax).';
                obj.y_arr = double(model.ymin:model.yres:model.ymax).';
                
                if length(input_fields) == 9
                    obj.z_arr = double(model.zmin:model.zres:model.zmax).';
                else
                    obj.z_arr = [];
                end
            case 6 % or the user may list them all out as arguments
                obj.x_arr = double(xmin:xres:xmax).';
                obj.y_arr = double(ymin:yres:ymax).';
                obj.z_arr = [];
            case 9 % if there is an imaging depth input
                obj.z_arr = double(zmin:zres:zmax).';
            otherwise
            error("Bad number of input arguments");
        end
    end
end
end