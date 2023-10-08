classdef DelayMatrix
% Delay Matrix class

    properties
        M           % a Delay matrix, in units of indices to delay.
        type        % they type of Delay Matrix
        Time        % how long it took to make the Delay Matrix
    end
    
    properties (Access = private)
        sens
        img
        data
        interpolate
        delays
        is_large = false;
    end

    properties (Dependent)
        M_gpu
    end
    
    properties (Access = private, Dependent)
        index_phase
    end

    methods
    
        function obj = DelayMatrix(sensor_array, image_area, dataset, version, interp)
        % obj = DelayMatrix(sensor_array, image_area, dataset, version, interp)
        %
        % Generates a Delay Matrix (one of several types) based on the Distance Formula for every
        % pixel location specified by image_area to every transducer
        % location specified in sensor_array.
        %
        %   Inputs:
        %       - sensor_array: A SensorArray object containing the
        %           locations of your transducers.
        %       - image_area: An ImageArea object containing the locations
        %           of the desired image pixels
        %       - dataset: A Dataset object containing rfdata and some
        %           other parameters
        %       - version (optional, default "index"): a string/character
        %           to specify which delay and sum schema this delay matrix
        %           is for. It can be "index", for DAS by indexing,
        %           "matrix" for DAS by sparse matrix multiplication, or
        %           "perspectiverotate" for delay matrix rotation DAS on a
        %           circular array.
        %       - interp (optional, default false): a bool true or false,
        %           specifying whether or not to interpolate between indices in
        %           the Delay Matrix. Iterpolation reduces computational
        %           time, but increases image fidelity slightly
        %
        %   Output:
        %       - obj.M: the delay matrix
        %       - obj.type: the type of delay matrix specified in 'version'

            arguments 
                sensor_array (1,1) SensorArray
                image_area (1,1) ImageArea
                dataset (1,1) Dataset
                version string {mustBeTextScalar, mustBeMember(version, ["index" "matrix" "perspectiverotate"])} = "index"
                interp (1,1) logical = false
            end

            tic;

            %Set class properties
            obj.type            = version;
            obj.interpolate   = interp;
            obj.sens          = sensor_array;
            obj.img           = image_area;
            obj.data          = dataset;
            
            %Run the subroutine that calculates delays based on sensor
            %array positions:
            obj                 = obj.distance_formula;
    
            switch version
                case "matrix"; obj = obj.delay_matrix_mmult;
                case "index"; obj = obj.delay_matrix_indexing;
                case "PerspectiveRotate"; obj = obj.delay_matrix_rotate;
                otherwise
                    error("Unrecognized version argument. Acceptable options are: 'matrix', 'index' (default), or 'PerspectiveRotate'");
            end
    
            obj.Time = toc;

        end
    
    end

    %subroutines
    methods (Access = private)
    
        function obj = distance_formula(obj)
            
            % Get the sensor array dimensions
            lxar                = obj.img.image_size(1);
            lyar                = obj.img.image_size(2);
            n_sens              = obj.sens.nSensors;
            
            % Argument checking for sizes
            ray_bytes   = ((lxar + lyar)*n_sens)*8;
            allr_bytes  = lxar*lyar*n_sens*2*8;
            gb          = 2^30;
    
            % If the matrix is larger than 20 GB, ask the user if they want
            % to continue
            if allr_bytes + ray_bytes > 20*gb
                warning("Requested RAM amount is large: "+ (allr_bytes + ray_bytes)/gb +" gB. Continue? (if yes, press any button, if not, press ctrl + c)")
                pause;
                obj.is_large = true;
            end
    
            % Calculate the delay matrix/indices
            xray(1,:,:)         = obj.img.x_arr*ones(1, n_sens)-ones(lxar,1)*obj.sens.x0;
            yray(:,1,:)         = obj.img.y_arr*ones(1, n_sens)-ones(lyar,1)*obj.sens.z0;
            obj.delays          = hypot(repmat(xray,lyar,1,1), repmat(yray,1,lxar,1)).*(1/obj.data.c/obj.data.dt);
    
        end
    
        function obj = delay_matrix_mmult(obj)
                            
            obj.type = "mmult";
            flrays   = cast(floor(obj.delays), 'uint32');
    
            % Get the row numbers (which simply repeats the sequence of 1:(lyar*lxar) lx0 times)
            I                   = flrays <= 0 & flrays > obj.data.frame_size;
            flrays(I)           = 1;
            Rows                = cast(repmat((1:obj.img.nPixels).', obj.sens.nSensors,1), 'single');
            Cols                = cast(flrays(:) + obj.index_phase(:), 'single');
    
            %Interpolation option. This adjusts the values at each index by
            %a factor of the truncated fraction of an index that would have otherwise
            %been rounded up or down. It adds the rest of that difference
            %to the next index.
            if obj.interpolate
                Rows            = [Rows;Rows];
                Cols            = [Cols;Cols+1];
                interprays      = obj.delays(:) - flrays(:);
                Vals            = cast([1-interprays; interprays],'double');
            else
                Vals            = ones(numel(obj.delays),1);
            end
    
            % Make the sparse matrix
            obj.M               = sparse(Rows, Cols, Vals, obj.img.nPixels, numel(obj.data.rfdata));
    
        end
    
        function obj = delay_matrix_indexing(obj)

            switch obj.data.pulse_ind
                case 0
                obj.delays = obj.delays;
    
                otherwise
                obj.delays = obj.delays + obj.data.pulse_ind;
            end
    
            obj.delays = cast(obj.delays,'uint32');
            obj.M = obj.delays + obj.index_phase;
            obj.M (obj.M > numel(obj.data.rfdata)) = numel(obj.data.rfdata);
            obj.M (obj.M < 1) = 1;

        end
        
        function obj = delay_matrix_rotate(obj)

            obj.type = "Index2D";

            % get size parameters
            frame_size          = uint32(obj.data.frame_size);
            lx0                 = length(obj.sens.x0);
            lxar                = obj.img.image_size(1);
            lyar                = obj.img.image_size(2);
            theta               = atan2(obj.sens.z0, obj.sens.x0);
            theta               = -theta*180/pi;

            % error check
            if lxar ~= lyar
                error("This method only supports square images.")
            end

            % find the imaging region dimensions
            Lx = abs(obj.img.x_arr(1)) + abs(obj.img.x_arr(end));
            Ly = abs(obj.img.y_arr(1)) + abs(obj.img.y_arr(end));
            dx = Lx/(lxar-1);
            dy = Ly/(lyar-1);
    
            % expand the imaging region so it's big enough to contain an
            % inscribed rotating circle of the center-to-corner line of the
            % square
            expansion_factor = sqrt(2);
            inflated.x_arr = obj.img.x_arr(1)*expansion_factor:dx:obj.img.x_arr(end)*expansion_factor;
            inflated.y_arr = obj.img.y_arr(1)*expansion_factor:dy:obj.img.y_arr(end)*expansion_factor;
            delay_matrix = uint32(sqrt((inflated.x_arr-obj.sens.x0(1)).^2+(inflated.y_arr.'-obj.sens.z0(1)).^2)/obj.data.c/obj.data.dt);

            % loop and rotate the image
            obj.M = zeros(lxar, lyar, lx0, 'uint32');
            for ii = 1:lx0
                phase = (ii-1)*frame_size;
                rotated_img = imrotate(delay_matrix, 90-theta(ii) ,'bilinear','crop');
                
                % index into the window of the actual imaging region size
                c_size = floor((size(rotated_img) - [lxar, lyar]) / 2);
                xind_vec = c_size(1):c_size(1) + lxar-1;
                yind_vec = c_size(2):c_size(2) + lyar-1;
                
                obj.M(:,:,ii) = phase + rotated_img(xind_vec, yind_vec);
            end

        end

    end

    methods

       function val = get.index_phase(delay_mat_obj)

            val(1,:,:)      = (0:size(delay_mat_obj.data.rfdata,2)-1)*(size(delay_mat_obj.data.rfdata,1));
            val             = repmat(val, length(delay_mat_obj.img.x_arr), length(delay_mat_obj.img.y_arr), 1);
            val             = cast(val, 'uint32');

       end

       function obj = set.sens(obj, val)

            if ~isa(val, "SensorArray")
                error("Expected a Sensor Array object as the sensor array input.");
            end

            obj.sens = val;
       end

       function obj = set.img(obj, val)
            
            if ~isa(val, "ImageArea")
                error("Expected an ImageArea object as the image region input");
            end
            obj.img = val;

       end

       function obj = set.data(obj, val)
            
            if ~isa(val, "Dataset")
                error("Expected a Dataset object for the data set input argument.");
            end
            obj.data = val;
       end

       function val = get.M_gpu(obj)

            val = gpuArray(obj.M);
    
       end

    end

end