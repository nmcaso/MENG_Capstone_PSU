classdef IndexMatrix

properties
    M
    Times
    R
    C
end

%primary constructor
methods

    function obj = IndexMatrix(sensorarray, imagearea, dataset, version, interp)

        switch version

           case "matrix"

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
%             Rows                = sort(Rows);
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
            obj.M               = sparse(Rows, Cols, Vals, lyar*lxar, numel(dataset.rfdata))./lx0;
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
            obj.Times.dims = toc;
            
            tic
            %2D distance matrices
            xray(1,:,:) = imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0;
            yray(:,1,:) = imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0;
            
            obj.Times.rays = toc;
            
            %3D index matrix
            allrays = hypot(repmat(xray,lyar,1,1), repmat(yray,1,lxar,1));
            
            switch dataset.pulse_ind
                case 0
                roundme = allrays.*(1/dataset.c/dataset.dt);

                otherwise
                roundme = allrays.*(1/dataset.c/dataset.dt) + dataset.pulse_ind;
            end

            allind           = cast(round(roundme),'uint32');
            obj.Times.allind = toc;
            
            tic
            %Add a constant to each page ofans allind for correct indexing later
            phasor1(1,:,:)  = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
            phasorcol       = repmat(phasor1,lxar,lyar,1);
            obj.Times.phasormat = toc;
            
            tic
            obj.M               = allind + cast(phasorcol, 'uint32');
            obj.Times.phaseadd  = toc;

            obj.Times.total = obj.Times.phaseadd+obj.Times.allind+obj.Times.rays+obj.Times.dims;
            
        end

    end

    function obj = Send2GPU(obj)
        obj.M   = gpuArray(obj.M);
    end

end

%alternate constructor
methods (Static)

    function obj = GPUIndexMatrix(sensorarray ,imagearea, dataset)
            imagearea.x_arr = gpuArray(imagearea.x_arr);
            imagearea.y_arr = gpuArray(imagearea.y_arr);
            sensorarray.x0  = gpuArray(sensorarray.x0);
            sensorarray.z0  = gpuArray(sensorarray.z0);

            tic
            %Sensor array dimensions
            lxar        = length(imagearea.x_arr);
            lx0         = length(sensorarray.x0);
            lyar        = length(imagearea.y_arr);
            lz0         = length(sensorarray.z0); 
            obj.Times.dims = toc;
            
            tic
            %2D distance matrices
            xray(1,:,:) = (imagearea.x_arr*ones(1, lx0,'gpuArray')-ones(lxar,1,'gpuArray')*sensorarray.x0).^2;
            yray(:,1,:) = (imagearea.y_arr*ones(1, lz0,'gpuArray')-ones(lyar,1,'gpuArray')*sensorarray.z0).^2;
            obj.Times.rays = toc;
            
            %3D index matrix
            allrays = repmat(xray,lyar,1,1) + repmat(yray,1,lxar,1);
            
            if dataset.pulse_ind == 0
                roundme = sqrt(allrays).*(1/dataset.c/dataset.dt);
            else
                roundme = sqrt(allrays).*(1/dataset.c/dataset.dt) + dataset.pulse_ind;
            end
            allind      = round(roundme);
            obj.Times.allind = toc;
            
            tic
            %Add a constant to each page of allind for correct indexing later
            smallphasor(1,:,:)  = (0:size(dataset.rfdata,2)-1)*(size(dataset.rfdata,1));
            phasor              = repmat(smallphasor,lxar,lyar,1);
            obj.Times.phasormat = toc;
            
            tic
            obj.M = allind+phasor;
            obj.Times.phaseadd = toc;

            obj.Times.total = obj.Times.phaseadd+obj.Times.allind+obj.Times.rays+obj.Times.dims;

    end

    function obj = ParIndexMatrix(sensorarray ,imagearea, dataset)
            
            tic
            %Sensor array dimensions
            lxar        = length(imagearea.x_arr);
            lx0         = length(sensorarray.x0);
            lyar        = length(imagearea.y_arr);
            lz0         = length(sensorarray.z0); 
            obj.Times.dims = toc;
            
            tic
            %2D distance matrices
            xray(1,:,:) = (imagearea.x_arr*ones(1, lx0)-ones(lxar,1)*sensorarray.x0).^2;
            yray(:,1,:) = (imagearea.y_arr*ones(1, lz0)-ones(lyar,1)*sensorarray.z0).^2;
            obj.Times.rays = toc;
            
            %3D index matrix
            allrays = repmat(xray,lyar,1,1) + repmat(yray,1,lxar,1);
            
            if dataset.pulse_ind == 0
                roundme = sqrt(allrays).*(1/dataset.c/dataset.dt);
            else
                roundme = sqrt(allrays).*(1/dataset.c/dataset.dt) + dataset.pulse_ind;
            end
            allind      = round(roundme);
            obj.Times.allind = toc;
            
            obj.M = allind;
            
            obj.Times.total = obj.Times.allind+obj.Times.rays+obj.Times.dims;

    end

end

end