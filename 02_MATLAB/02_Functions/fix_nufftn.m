function y = fix_nufftn(x, t, f)
% NUFFTN  N-dimensional nonuniform Discrete Fourier Transform.
%
%   y = nufftn(x,t) approximates the N-dimensional nonuniform discrete
%   Fourier transform (NUDFT) of x using the sample points specified by t.
%   The sample points t can be specified as an cell array of grid vectors
%   for each dimension, or as an N by D coordinate matrix, where t(i,:) is
%   a set of D-dimensional coordinates corresponding to the input value
%   x(i).
%
%   The query points of the transform are a uniform D-dimensional grid of
%   M = ceil(numel(x)^(1/D)) points in each dimension  The output y is a
%   vector of M^d transform coefficients corresponding to the transform
%   values at each query point.
%
%   y = nufftn(x,t,f) approximates the N-dimensional NUDFT of x using the
%   sample points specified by t and the query points specified by f.  The
%   query points f can be specified as an cell array of grid vectors for
%   each dimension, or as an M by D coordinate matrix, where f(i,:) is a
%   set of D-dimensional coordinates corresponding to the output value y(i).
%
%   The output y is a length-M vector of transform coefficients, where M is
%   the total number of query points.
%
%   Example: Compute a 2-D DFT of uniform data onto a set of spiral points.
%       x = peaks(128);
%       % Create a set of spiral points for the transform.
%       ang = linspace(0, 10, 1024)';
%       mag = linspace(0, 0.5, 1024)';
%       f = mag .* [cospi(ang), sinpi(ang)];
%       % Compute the uniform-to-non-uniform DFT.
%       X = nufftn(x, [], f);
%       % Visualize the magnitude of the results.
%       scatter3(f(:,1), f(:,2), log10(abs(X)));
%
%   Example: Compute an 3-D DFT of non-uniform data onto a set of points.
%       % Non-uniform source.
%       t = (rand(4096, 3) * 16) - 8;
%       x = exp(-(t(:,1).^2 + t(:,2).^2 + t(:,3).^2)/16);
%       % Uniform nodes.
%       f = { (0:15)/16, (0:15)/16 (0:15)/16 };
%       X = nufftn(x, t, f);
%
%   See also nufft fftn ifftn

%  Copyright 2019-2022 The MathWorks, Inc.

    % Parse the input array.
    if isinteger(x)
        if isa(x, 'uint64') || isa(x, 'int64')
            error(message('MATLAB:fftfcn:InvalidInputType'));
        end
        x = double(x);
    elseif islogical(x)
        x = double(x);
    elseif ~isfloat(x)
        error(message('MATLAB:fftfcn:InvalidInputType'));
    elseif issparse(x)
        error(message('MATLAB:nufft:nufftNoSparseArrays'));
    end
    
    % Set default sample and query points if not specified.
    if nargin < 2
        t = [];
    end
    if nargin < 3
        f = [];
    end
    
    % Default nodes collapse to FFTN.
    if isempty(t) && isempty(f)
        y = fftn(x);
        return;
    end
    
    % Determine the number of dimensions.
    dims = [ndims(x), 0, 0];
    if ~isempty(t)
        % Determine the dimensionality of the sample points.
        if iscell(t)
            dims(2) = numel(t);
        else
            dims(2) = size(t, 2);
        end
    end
    if ~isempty(f)
        % Determine the dimensionality of the query points.
        if iscell(f)
            dims(3) = numel(f);
        else
            dims(3) = size(f, 2);
        end
    end
    % The number of dimensions is the maximum between all of the prescribed
    % dimensions.
    dims = max(dims);
    
    % Parse the node inputs.
    isSamplePoints = true;
    snodes = parseNodes(t, dims, x, isSamplePoints);
    isSamplePoints = false;
    qnodes = parseNodes(f, dims, x, isSamplePoints);
    
    % Check that the sample points match.
    if prod(snodes.count) ~= numel(x)
        if snodes.type == "gridded"
            error(message('MATLAB:nufft:nufftnSamplePointsMismatch'));
        else
            error(message('MATLAB:nufft:nufftnScatteredSamplePointsMismatch'));
        end
    end
    
    % Check for simple cases.
    if (snodes.type ~= "scattered") && (qnodes.type ~= "scattered")
        % Separable NUFFT.
        y = nufftn_separable(x, snodes, qnodes);
    elseif (~snodes.uniform && ~qnodes.uniform)
        % Type-IV NUFFT.
        y = matlab.internal.math.nufftndirect(x, snodes.values, ...
            qnodes.values);
    else
        % Type-II or Type-III NUFFT.
        y = nufftn_interpolate(x, snodes, qnodes);
    end
    % Always return a column.
    y = y(:);
    
end

%--------------------------------------------------------------------------
function nodeData = parseNodes(pts, dims, x, isSamplePoints)
% Parse node input data.
    % Default nodes.
    nodeData.flip = false(1, dims);
    if isequal(pts, [])
        % Implicit nodes match the size of x.
        nodeData.type = "implicit";
        nodeData.count = ones(1, dims);
        nodeData.uniform = true;
        if isSamplePoints
            nodeData.count(1:ndims(x)) = size(x);
        else
            nodeData.count = nodeData.count * ceil(numel(x).^(1/dims));
        end
        nodeData.values = [];
    elseif iscell(pts)
        % Gridded data.
        nodeData.type = "gridded";
        nodeData.values = repmat({ zeros(1,1,class(x)) }, 1, dims);
        nodeData.count = ones(1, dims);
        nodeData.uniform = true;
        cls = class(x);
        for i = 1:length(pts)
            % Validate each node grid.
            if ~isfloat(pts{i})
                error(message('MATLAB:nufft:invalidGridVectorDataType'));
            elseif ~isvector(pts{i}) || ...
                isempty(pts{i}) || ~allfinite(pts{i}) || ...
                ~isreal(pts{i}) || issparse(pts{i})
                error(message('MATLAB:nufft:ndGridVectorsNotFullReal'));
            end
            tf = matlab.internal.math.isuniform(pts{i});
            if ~tf && matlab.internal.math.isuniform(flip(pts{i}))
                tf = true;
                pts{i} = flip(pts{i});
                nodeData.flip(i) = true;
            end
            nodeData.uniform = nodeData.uniform && tf;
            nodeData.values{i} = cast(pts{i}, cls);
            nodeData.count(i) = numel(pts{i});
        end
    else
        % Matrix data.
        if ~isfloat(pts)
            error(message('MATLAB:nufft:invalidNodeDataType'));
        elseif ~ismatrix(pts) || ...
           ~allfinite(pts) || ~isreal(pts) || issparse(pts)
            error(message('MATLAB:nufft:ndScatteredNodesNotFullReal'));
        end
        nodeData.count = size(pts, 1);
        if nodeData.count == 1
            % Treat single points as gridded data.
            nodeData.count = ones(1, dims);
            nodeData.type = "gridded";
            nodeData.values = repmat({ zeros(1,1,class(x)) }, 1, dims);
            for dim = 1:size(pts,2)
                nodeData.values{dim} = cast(pts(dim), class(x));
            end
            nodeData.uniform = true;
        else
            nodeData.type = "scattered";
            nodeData.uniform = false;
            nodeData.values = [pts, ...
                zeros(nodeData.count, dims-size(pts,2), class(x))];
        end
    end
end

%--------------------------------------------------------------------------
function y = nufftn_separable(x, snodes, qnodes)
% Perform the N-D NUFFT on a separable set of grids.

    % Make sure x is reshaped to the correct size.
    if snodes.type == "gridded"
        sz = cellfun(@length, snodes.values, 'UniformOutput', true);
        x = reshape(x, sz);
    end
    % Optimize the dimensional order.
    [~,dimorder] = sort(qnodes.count ./ snodes.count);
    y = x;
    for d = 1:length(dimorder)
        % Get the sample points for this dimension.
        t = [];
        if snodes.type == "gridded"
            t = snodes.values{dimorder(d)};
            if snodes.flip(dimorder(d))
                y = flip(y, dimorder(d));
            end
        end
        % Get the query points for this dimension.
        if qnodes.type == "implicit"
            m = qnodes.count(dimorder(d));
            f = (0:(m-1))/m;
        else
            f = qnodes.values{dimorder(d)};
        end
        % Invoke the 1-D NUFFT on the given dimension.
        y = nufft(y, t, f, dimorder(d));
        if qnodes.flip(dimorder(d))
            y = flip(y, dimorder(d));
        end
    end
end

%--------------------------------------------------------------------------
function y = nufftn_interpolate(x, snodes, qnodes)
% N-D NUFFT via interpolation.
    dims = max(numel(qnodes.count), numel(snodes.count));
    if qnodes.type == "scattered"
        % Type-II NUFFTN.
        % Query nodes will lie in [0 1)
        % Sample points will lie in [0 N)
        needsShift = snodes.type == "gridded";
        if needsShift
            % Sample points are uniformly gridded.
            needsShift = true;
            forig = qnodes.values;
            [qnodes, N, phaseOffsets] = shiftAndScaleNodes(qnodes, snodes);
            if any(snodes.flip)
                % Apply flips to the input nodes.
                x = reshape(x, snodes.count);
                for d = find(snodes.flip)
                    x = flip(x, d);
                end
            end
        else
            qnodes.values = mod(qnodes.values, 1);
            N = snodes.count;
        end
        % Squeeze out singleton dimensions.
        nonSingletonDims = N > 1;
        if any(nonSingletonDims)
            N = N(nonSingletonDims);
            qnodes.values = qnodes.values(:,nonSingletonDims);
            % Select interpolation parameters.
            [n, gparams] = selectInterpolationParameters(N);
            % Compute the transform.
            y = matlab.internal.math.nufftninterp(reshape(x, [N, 1]), n, ...
                [], qnodes.values, gparams);
        else
            y = repmat(x, size(qnodes.values, 1), 1);
        end
        % Apply the phase shift.
        if needsShift
            for dim = 1:dims
                y = matlab.internal.math.unitPhaseFactor(...
                    forig(:,dim), phaseOffsets(dim)) .* y;
            end
        end
    else
        % Type-III NUFFTN.
        % Query nodes will lie in [0 N)
        % Sample points will lie in [0 1)
        needsShift = qnodes.type == "gridded";
        if needsShift
            % Gridded nodes, but uniform.
            torig = snodes.values;
            [snodes, N, phaseOffsets] = shiftAndScaleNodes(snodes, qnodes);
            % Apply the pre-scaling.
            for dim = 1:dims
                x = matlab.internal.math.unitPhaseFactor(...
                    torig(:,dim), phaseOffsets(dim)) .* x(:);
            end
        else
            N = qnodes.count;
            snodes.values = mod(snodes.values, 1); 
        end
        % Squeeze out singleton dimensions.
        nonSingletonDims = N > 1;
        if any(nonSingletonDims)
            N = N(nonSingletonDims);
            snodes.values = snodes.values(:,nonSingletonDims);
            % Select interpolation parameters.
            [n, gparams] = selectInterpolationParameters(N);
            % Compute the transform.
            y = matlab.internal.math.nufftninterp(x, n, ...
                snodes.values, N, gparams);
        else
            y = sum(x(:));
        end
        if needsShift && any(qnodes.flip)
            y = reshape(y, qnodes.count);
            for d = find(qnodes.flip)
                y = flip(y, d);
            end
            y = y(:);
        end
    end
end

%--------------------------------------------------------------------------
function [nuNodes, N, phaseOffs] = shiftAndScaleNodes(nuNodes, gridNodes)
% Compute phase shift and scale non-uniform nodes to turn uniform gridded
% nodes into implicit nodes.
    dims = max(numel(nuNodes.count), numel(gridNodes.count));
    % Compute the phase offsets, if there are any.
    phaseOffs = zeros(dims, 1);
    % Gridded nodes, but uniform.
    N = ones(1, dims);
    for dim = 1:dims
        % Determine the grid size.
        N(dim) = numel(gridNodes.values{dim});
        % Shift the grid to be uniform.
        phaseOffs(dim) = min(gridNodes.values{dim});
        if N(dim) > 1
            delta = median(diff(gridNodes.values{dim}));
        else
            % No need to rescale scalar data.
            delta = 1;
        end
        % Rescale the non-uniform nodes by the uniform grid spacing, then
        % apply periodicity.
        nuNodes.values(:,dim) = mod(nuNodes.values(:,dim) * delta, 1);
    end
end

%--------------------------------------------------------------------------
function [n, gparams] = selectInterpolationParameters(N)
% Select the interpolation parameters for a type-II or type-III N-D NUFFT.
    N = N(:);
    % Interpolatory window length grows logarithmically with input size.
    m = 4 + log(N);
    n = ceil(1.75.*N);
    % Dimensions whose sizes are too small should be handled directly.
    smallDims = N < 16;
    n(smallDims) = N(smallDims);
    gparams = [max(2,m/2.5), m./n];
    gparams(smallDims,2) = 0.5;
    n = n';
end