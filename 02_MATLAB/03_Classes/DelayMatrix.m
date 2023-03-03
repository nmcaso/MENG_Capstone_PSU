classdef DelayMatrix
% Delay Matrix class

properties
    M           % the delay matrix
    Times       % for profiling
    type        % an indicator for some of the DAS algorithms
end

methods

    function delay_mat_obj = DelayMatrix(sensorarray, imagearea, dataset, version, interp)
    % constructor - generates a delay matrix of several types

        switch version
        % valid versions are "matrix", "index", and "PerspectiveRotate"

            case "matrix"
            % generates a sparse matrix for the matrix mult method

            delay_mat_obj.type = "mmult";
            if ~isempty(imagearea.z_arr)
                warning("Matrix Multiplication unsupported for depth arguments. Only performing image reconstruction in 2D")
            end

            % Get the sensor array dimensions
            tic
            lxar                = length(imagearea.x_arr);
            lx0                 = length(sensorarray.x0);
            lyar                = length(imagearea.y_arr);
            lz0                 = length(sensorarray.z0);
            delay_mat_obj.Times.dims      = toc;

            % Calculate the delay matrix/indices
            tic
            xray(1,:,:)         = imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0;
            yray(:,1,:)         = imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0;
            allrays             = hypot(repmat(xray,lyar,1,1), repmat(yray,1,lxar,1)).*(1/dataset.c/dataset.dt);
            flrays              = floor(allrays);
            delay_mat_obj.Times.rays      = toc;

            % Get the row numbers (which simply repeats the sequence of 1:(lyar*lxar) lx0 times)
            tic
            I                   = flrays > 0 & flrays < size(dataset.rfdata,1);
            flrays(~I)          = 1;
            Rows                = cast(repmat((1:lxar*lyar).',lx0,1),'single');
            delay_mat_obj.Times.getrows   = toc;

            % Get the column numbers from the indices
            tic
            phasor1(1,:,:)      = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
            phasorcol           = repmat(phasor1,lxar,lyar,1);
            Cols                = cast(flrays(:) + phasorcol(:),'single');
            delay_mat_obj.Times.getcols   = toc;

            %Interpolation option. This adjusts the values at each index by
            %a factor of the truncated fraction of an index that would have otherwise
            %been rounded up or down. It adds the rest of that difference
            %to the next index.
            if interp
                Rows            = [Rows;Rows];
                Cols            = [Cols;Cols+1];
                interprays      = allrays(:) - flrays(:);
                Vals            = cast([1-interprays; interprays],'double');
            else
                Vals            = ones(numel(allrays),1);
            end

            % Make the sparse matrix
            tic
            delay_mat_obj.M               = sparse(Rows, Cols, Vals, lyar*lxar, numel(dataset.rfdata));
            delay_mat_obj.Times.matrix    = toc;

            % Capture Performance data
            delay_mat_obj.Times.total     = delay_mat_obj.Times.matrix+delay_mat_obj.Times.getcols+delay_mat_obj.Times.getrows+delay_mat_obj.Times.rays;

            case "index"
            % This option generates a 3D index matrix

            tic
            % Get the sensor array dimensions
            lxar        = length(imagearea.x_arr);
            lx0         = length(sensorarray.x0);
            lyar        = length(imagearea.y_arr);
            lz0         = length(sensorarray.z0);
            lzar        = length(imagearea.z_arr);

            % Argument checking for sizes
            ray_bytes = (lxar*lx0 + lyar*lz0)*8;
            allr_bytes = lxar*lyar*lx0*2*8;
            bigflag = false;

            % If the matrix is larger than 20 GB, ask the user if they want
            % to continue
            if allr_bytes + ray_bytes > 20*2^30
                warning("Requested RAM amount is large: "+ (allr_bytes + ray_bytes)/2^30 +" gB. Continue? (if yes, press any button, if not, press ctrl + c)")
                pause;
                bigflag = true;
            end
            delay_mat_obj.Times.dims = toc;
            
            %2D distance matrices
            tic
            if lzar == 0
                xray(1,:,:) = imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0;
                yray(:,1,:) = imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0;
                delay_mat_obj.type = "Index2D";
            else
                xray(1,:,:,:) = repmat(imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0, [1,1,1,lzar]);
                yray(:,1,:,:) = repmat(imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0, [1,1,1,lzar]);
                zray(1,1,:,:) = (imagearea.z_arr*ones(1, lx0)-ones(lzar,1)*zeros(1, lz0)).';
                
                xray = repmat(xray,[lyar,1,1,1]);
                yray = repmat(yray,[1,lxar,1,1]);
                zray = repmat(zray,[lxar,lyar,1,1]);
                delay_mat_obj.type = "Index3D";
            end
            if bigflag
                xray = single(xray);
                yray = single(yray);
            end
            delay_mat_obj.Times.rays = toc;
            
            %3D delay matrix
            if lzar == 0
                allrays = hypot(repmat(xray,lyar,1,1), repmat(yray,1,lxar,1));
            else
                allrays = zeros([size(xray,[1 2]) size(xray,4) size(xray, 3)]);
            for ii = 1:lx0
                for jj = 1:lzar
                    allrays(:,:,jj,ii) = sqrt(xray(:,:,ii, jj).^2 + yray(:,:,ii,jj).^2) + zray(:,:,ii,jj); 
                end
            end
            end

            switch dataset.pulse_ind
                case 0
                roundme = allrays.*(1/dataset.c/dataset.dt);

                otherwise
                roundme = allrays.*(1/dataset.c/dataset.dt) + dataset.pulse_ind;
            end

            allind           = cast(roundme,'uint32');
            delay_mat_obj.Times.allind = toc;
            
            %Add a constant to each page of of the delay matrix for correct indexing later
            tic
            if lzar == 0
                phasor1(1,:,:)  = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
                phasorcol       = repmat(phasor1,lxar,lyar,1);
            else
                phasor1(1,1,:,:)  = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
                phasorcol       = repmat(phasor1,[lxar,lyar,lzar,1]);
            end            
            delay_mat_obj.Times.phasormat = toc;
            
            tic
            delay_mat_obj.M               = allind + cast(phasorcol, 'uint32');
            delay_mat_obj.M(delay_mat_obj.M > numel(dataset.rfdata)) = numel(dataset.rfdata);
            delay_mat_obj.M(delay_mat_obj.M < 1) = 1;
            delay_mat_obj.Times.phaseadd  = toc;

            % Output the time performance data
            delay_mat_obj.Times.total = delay_mat_obj.Times.phaseadd+delay_mat_obj.Times.allind+delay_mat_obj.Times.rays+delay_mat_obj.Times.dims;
            
            case "PerspectiveRotate"
            %generates a perspective-rotated 3D index matrix

            delay_mat_obj.type = "Index2D";
            tic;

            % get size parameters
            frame_size          = uint32(size(dataset.rfdata,1));
            lxar                = length(imagearea.x_arr);
            lx0                 = length(sensorarray.x0);
            lyar                = length(imagearea.y_arr);
            [theta, ~]          = cart2pol(sensorarray.x0, sensorarray.z0);
            theta               = -theta*180/pi;

            % error check
            if lxar ~= lyar
                error("This method is only supported for square arrays currently")
            end

            % find the imaging region dimensions
            Lx = abs(imagearea.x_arr(1)) + abs(imagearea.x_arr(end));
            Ly = abs(imagearea.y_arr(1)) + abs(imagearea.y_arr(end));
            dx = Lx/(lxar-1);
            dy = Ly/(lyar-1);
    
            % expand the imaging region so it's big enough to contain an
            % inscribed rotating circle of the center-to-corner line of the
            % square
            expansion_factor = sqrt(2);
            inflated.x_arr = imagearea.x_arr(1)*expansion_factor:dx:imagearea.x_arr(end)*expansion_factor;
            inflated.y_arr = imagearea.y_arr(1)*expansion_factor:dy:imagearea.y_arr(end)*expansion_factor;
            delay_matrix = uint32(sqrt((inflated.x_arr-sensorarray.x0(1)).^2+(inflated.y_arr.'-sensorarray.z0(1)).^2)/dataset.c/dataset.dt);

            % loop and rotate the image
            delay_mat_obj.M = zeros(lxar, lyar, lx0, 'uint32');
            for ii = 1:lx0
            phase = (ii-1)*frame_size;
            rotated_img = imrotate(delay_matrix, 90-theta(ii) ,'bilinear','crop');
            
            % index into the window of the actual imaging region size
            c_size = floor((size(rotated_img) - [lxar, lyar]) / 2);
            xind_vec = c_size(1):c_size(1) + lxar-1;
            yind_vec = c_size(2):c_size(2) + lyar-1;
            
            delay_mat_obj.M(:,:,ii) = phase + rotated_img(xind_vec, yind_vec);
            end

            delay_mat_obj.Times.total = toc;
        end

    end

    function obj = Send2GPU(obj)
    %places the delay matrix onto the GPU
        
        obj.M   = gpuArray(obj.M);
        obj.type = "gpu" + obj.type;
    end

end

end