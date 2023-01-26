%Imaging Region

classdef ImageArea

    properties
        x_arr
        y_arr
        z_arr
        pulse_ind
    end
    
    methods
    
        %Constructor
        function obj = ImageArea(xmax, xmin, xres, ymax, ymin, yres, zmax, zmin, zres)
            
            switch nargin

                case 1
                    model = xmax;
                    if ~isa(model,'struct')
                        error("The single input option for the class constructor requires a structure.")
                    end

                    input_fields = string(fields(orderfields(model)));
                    predefined_fields = ["xmax"; "xmin"; "xres"; "ymax"; "ymin"; "yres"; "zmax"; "zmin"; "zres"];
                    
                    if mod(length(input_fields),3) ~= 0 || ~all(input_fields == predefined_fields(1:length(input_fields)))
                        error("Bad input arguments. The structure fields must be: xmin, xres, xmax, ymin, yres, ymax, and optionally zmin, zres, zmax");
                    end

                    obj.x_arr = (model.xmin:model.xres:model.xmax).';
                    obj.y_arr = (model.ymin:model.yres:model.ymax).';
                    
                    if length(input_fields) == 9
                        obj.z_arr = (model.zmin:model.zres:model.zmax).';
                    else
                        obj.z_arr = [];
                    end

                case 6
                
                obj.x_arr = (xmin:xres:xmax).';
                obj.y_arr = (ymin:yres:ymax).';
                obj.z_arr = [];

                case 9
                    
                obj.z_arr = (zmin:zres:zmax).';

                otherwise
                    
                error("Bad number of input arguments");

            end
    
        end

    end

end