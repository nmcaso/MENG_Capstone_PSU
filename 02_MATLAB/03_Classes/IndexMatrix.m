classdef IndexMatrix

properties
    M
    Times
    R
    C
    type
end

%primary constructor
methods

    function obj = IndexMatrix(sensorarray, imagearea, dataset, version, interp)

        switch version

            case "matrix"
            obj.type = "mmult";
            if ~isempty(imagearea.z_arr)
                warning("Matrix Multiplication unsupported for depth arguments. Only performing image reconstruction in 2D")
            end

            %Sensor array dimensions
            tic
            lxar                = length(imagearea.x_arr);
            lx0                 = length(sensorarray.x0);
            lyar                = length(imagearea.y_arr);
            lz0                 = length(sensorarray.z0);
            obj.Times.dims      = toc;

            %Calculate the delay matrix/indices
            tic
            xray(1,:,:)         = imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0;
            yray(:,1,:)         = imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0;
            allrays             = hypot(repmat(xray,lyar,1,1), repmat(yray,1,lxar,1)).*(1/dataset.c/dataset.dt);
            flrays              = floor(allrays);
            obj.Times.rays      = toc;

            %Get the row numbers (which simply repeats the sequence of 1:(lyar*lxar) lx0 times)
            tic
            I                   = flrays > 0 & flrays < size(dataset.rfdata,1);
            flrays(~I)          = 1;
            Rows                = cast(repmat((1:lxar*lyar).',lx0,1),'single');
            obj.Times.getrows   = toc;

            %Get the column numbers from the indices
            tic
            phasor1(1,:,:)      = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
            phasorcol           = repmat(phasor1,lxar,lyar,1);
            Cols                = flrays(:) + phasorcol(:);
            obj.Times.getcols   = toc;

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

            %Generate the sparse DAS matrix
            tic
            obj.M               = sparse(Rows, Cols, Vals, lyar*lxar, numel(dataset.rfdata));
            obj.Times.matrix    = toc;

            %Capture Performance data
            obj.Times.total     = obj.Times.matrix+obj.Times.getcols+obj.Times.getrows+obj.Times.rays;

            case "index"

            tic
            %Sensor array dimensions
            lxar        = length(imagearea.x_arr);
            lx0         = length(sensorarray.x0);
            lyar        = length(imagearea.y_arr);
            lz0         = length(sensorarray.z0);
            lzar        = length(imagearea.z_arr);

            %Argument checking
            ray_bytes = (lxar*lx0 + lyar*lz0)*8;
            allr_bytes = lxar*lyar*lx0*2*8;
            bigflag = false;

            if allr_bytes + ray_bytes > 10*2^30
                warning("Requested RAM amount is large: "+ (allr_bytes + ray_bytes)/2^30 +" gB. Continue?")
                pause;
                bigflag = true;
            end

            obj.Times.dims = toc;
            
            tic
            %2D distance matrices
            if lzar == 0
                xray(1,:,:) = imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0;
                yray(:,1,:) = imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0;
                obj.type = "Index2D";
            else
                xray(1,:,:,:) = repmat(imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0, [1,1,1,lzar]);
                yray(:,1,:,:) = repmat(imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0, [1,1,1,lzar]);
                zray(1,1,:,:) = (imagearea.z_arr*ones(1, lx0)-ones(lzar,1)*zeros(1, lz0)).';
                
                xray = repmat(xray,[lyar,1,1,1]);
                yray = repmat(yray,[1,lxar,1,1]);
                zray = repmat(zray,[lxar,lyar,1,1]);
                obj.type = "Index3D";
            end
            if bigflag
                xray = single(xray);
                yray = single(yray);
            end

            obj.Times.rays = toc;
            
            %3D index matrix

            if lzar == 0
            allrays = hypot(repmat(xray,lyar,1,1), repmat(yray,1,lxar,1));
            else
            allrays = sqrt(xray.^2 + yray.^2 + zray.^2); % questionably correct
            allrays = permute(allrays, [1 2 4 3]);
            end

            switch dataset.pulse_ind
                case 0
                roundme = allrays.*(1/dataset.c/dataset.dt);

                otherwise
                roundme = allrays.*(1/dataset.c/dataset.dt) + dataset.pulse_ind;
            end

            allind           = cast(roundme,'uint32');
            obj.Times.allind = toc;
            
            tic
            %Add a constant to each page ofans allind for correct indexing later
            if lzar == 0
            phasor1(1,:,:)  = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
            phasorcol       = repmat(phasor1,lxar,lyar,1);
            else
            phasor1(1,1,:,:)  = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
            phasorcol       = repmat(phasor1,[lxar,lyar,lzar,1]);
            end            
            obj.Times.phasormat = toc;
            
            tic
            obj.M               = allind + cast(phasorcol, 'uint32');
            obj.M(obj.M > numel(dataset.rfdata)) = numel(dataset.rfdata);
            obj.Times.phaseadd  = toc;

            obj.Times.total = obj.Times.phaseadd+obj.Times.allind+obj.Times.rays+obj.Times.dims;
            
        end

    end

    function obj = Send2GPU(obj)
        obj.M   = gpuArray(obj.M);
        obj.type = "gpu" + obj.type;
    end


end

end