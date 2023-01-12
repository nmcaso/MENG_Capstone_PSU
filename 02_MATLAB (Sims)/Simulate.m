            %% Setup grid and configuration structures
            clear; clc; tau = 2*pi; close all;                              %tau > pi
workdir     = 'C:\Users\cason\OneDrive\Documents\PSU\Project\02_MATLAB (Sims)\'; %working directory
refdir      = [workdir(1:55) '\03_Classes\'];                               %import DAS code classes
cd          (workdir);                                                      
addpath     (refdir);
load        ("01_DataIn\DatasetSample.mat");                                %Load sample sensor configurations
load        ("01_DataIn\SensorArraySample.mat");                            
% solvetype   = 1;                                                            %finite-difference in 3D
solvetype   = 2;                                                            %K-Space method
% [sensarr.x0, sensarr.z0] = sensordata();

            %load configuration variables
config.     dt          = 1/(19.2308*1E6); % source:leafdata.fs             % [s] transducer time-spacing
config.     N           = 300; %size(leafdata.rfdata,1);                          % number of points per DAS imaging 'frame'
config.     T           = config.N*config.dt;                               % [s] total time for 1 frame
config.     t           = 0:config.dt:config.T;                             % [s] time vector
config.     f           = 1e5;                                              % [Hz] pulsed laser frequency
config.     omg         = tau*config.f;                                     % [rad] pulsed laser angular frequency
config.     c           = 1.4895E3; %source: leafdata.c                     % [m/s] inhomogeneous speed of sound
config.     lamda       = config.c/config.f;                                % [m] pulsed laser wavelength
config.     pulse_len   = 1200;                                             % number of samples to pulse the laser for
config.     rho         = 997;                                              % [kg/m^2] tissue density
config.     animate     = true;                                            % switch to plot graphs over time
config.     complexp    = false;                                            % switch to use complex pressure or real pressure

            %model boundary layer switches
model.      zboundaries = true;                                             % switch for z boundaries (false = periodic)
model.      PML         = true;                                             % switch for absorbing layer boundaries (false = periodic)

            %model dimensions
model.      widthfac    = 4;                                                % (model width)/(size of transducer)
model.      dens        = 15;                                                % number of points per wavelength in all directions
model.      depth       = 12e-3;                                            % [m] mean depth of upper arm skin thickness for male with BMI 24
model.      length      = model.widthfac*(max(sensarr.x0)-min(sensarr.x0)); % [m] width (and length)
model.      resolutions = [ config.lamda/(model.dens);                      % [m] dx
                            config.lamda/(model.dens);                      % [m] dy
                            config.lamda/(model.dens)  ];                   % [m] dz
model.      size        = [ ceil(model.length/model.resolutions(1)), ...    % [m] L_x
                            ceil(model.length/model.resolutions(2)), ...    % [m] L_y
                            ceil(model.depth/model.resolutions(3))];        % [m] L_z
model.      cfl         = config.c*config.dt./model.resolutions;            % CFL number

            %absorbing layer properties
            if model.PML
model.      Pml.Depth   = 40;                                               % number of points for the absorbing layer
model.      Pml.Window  = hann(2*model.Pml.Depth,'periodic');               % multiplier function for the window
model.      Pml.halfx   = model.Pml.Window(1:model.Pml.Depth);              % window in the x-direction
model.      Pml.halfy(1,:,:) = model.Pml.Window(1:model.Pml.Depth);         % window in the y-direction
            end
            
disp("Model Size" + [" x: " " y: " " z: "] + string(model.size));           % print model sizes to the console
disp("CFL" + [" x "," y "," z "].' + model.cfl)                             % print CFL to the console

            %% Set up a Source and time vectors
source      .strength   = 0.1;                                              % [Pa] non-complex amplitude
            switch config.complexp
                case true
pfun        = @(t) source.strength*exp(1j*config.omg*t);                    % complex pressure function for point source
                case false
pfun        = @(t) source.strength*cos(config.omg*t);                       % non-complex pressure function for point source
            end

            % model a blood vessel with unitary radius
source.     X = (0.8*model.length:-model.resolutions(1):0.1*model.length).';% [m] x coordinates for source
source.     Y = cos(source.X./max(source.X))*0.8*model.length;              % [m] y coordinates for source
source.     Z = ones(length(source.Y),1)*model.depth/2;                     % [m] z coordinates for source

            %model a point-source
% source.     X           = model.length/2;
% source.     Y           = model.length/2;
% source.     Z           = model.depth/2;

            % model a spherical "blob" source
% blobradius  = 0.005;
% blobcenterxy= model.length/2;
% blobcenterz = model.depth/2;
% [sx, sy, sz] = sphere(50);
% source.     X           = 0.5*(1+sx(:)).*blobradius + blobcenterxy;
% source.     Y           = 0.5*(1+sy(:)).*blobradius + blobcenterxy;
% source.     Z           = 0.5*(1+sz(:)).*blobradius + blobcenterz;

            %calculate source indices
source.     indx        = round(source.X/model.resolutions(1));
source.     indy        = round(source.Y/model.resolutions(2));
source.     indz        = round(source.Z/model.resolutions(3));
% source.     indz        = ones(size(source.indy))*2;
source.     ind         = sub2ind(model.size,source.indx, ...
                        source.indy, source.indz);                          %calculate linear indices for the model matrix

            %Set 2 pressure fields for initial conditions
source.     p0          = zeros(model.size);
source.     p0(source.ind) = pfun(config.t(1));                             %[Pa] pressure 0 is the source strength as a monopole at time 0
source.     p1          = source.p0;                                        %[Pa] pressure 1 equal to pressure 0

            % set the pressure at t = 1
            for wv = 1:2                                                    % for the forward and backward-going waves
            %get the indices of the adjacent points to the source
xsubs       = sub2ind(model.size,source.indx+2*wv-3,source.indy,source.indz);
ysubs       = sub2ind(model.size,source.indx,source.indy+2*wv-3,source.indz);
zsubs       = sub2ind(model.size,source.indx,source.indy,source.indz+2*wv-3);

            % propagate the pressure 1 time step forward for the adj. points.
source.     p1(xsubs)   = model.cfl(1)*pfun(config.t(2));                   %[Pa]
source.     p1(ysubs)   = model.cfl(2)*pfun(config.t(2));                   %[Pa]
source.     p1(zsubs)   = model.cfl(3)*pfun(config.t(2));                   %[Pa]a
            end
            
            % and increment the pressure at the point source(s)
source.     p1(source.ind) = pfun(config.t(2));                             %[Pa]

            %set the speed of sound-contrast to match boundary conditions
            %(for k-space method)
source.     cxyz                    = ones(model.size).*config.c;           % 1498 [m/s] tissue speed of sound
source.     cxyz(:,:,1)             = 343;                                  % [m/s] air speed of sound at top of model
source.     cxyz(:,:,end)           = 2105;                                 % [m/s] bone speed of sound at bottom of model

            %set the density contrast to match b.c's
source.     rhoxyz                  = ones(model.size).*config.rho;         % [kg/m^3]
source.     rhoxyz(:,:,1)           = 1.293;                                % [kg/m^3] air density
source.     rhoxyz(:,:,end)         = 1900;                                 % [kg/m^3] bone density
source.     laplacianrho            = zeros(model.size);                   
source.     laplacianrho(:,:,1)     = 1; %(1.293 - config.rho);             % a guess at the laplacian of density
source.     laplacianrho(:,:,end)   = 1; %(1900 - config.rho);

            %% Define observation points
probe.      actualx     = sensarr.x0 + (model.length)/2;                    % [m] x coordinates of actual transducers
probe.      actualy     = sensarr.z0 + (model.length)/2;                    % [m] y coordinates of transducers

            %calculate indices from the coordinates
probe.      indx        = round(probe.actualx/model.resolutions(1));
probe.      indy        = round(probe.actualy/model.resolutions(2));
probe.      indz        = ones(1,length(probe.indy))*(model.size(3) - 1);   % z = 1 is the bottom of the model.

            %catch indices that are beyond the model size boundaries
flds        = fields(probe);
            for dim = flds(3:5)'
probe.(dim{1})(probe.(dim{1}) == 0) = 1;
probe.(dim{1})(probe.(dim{1}) > model.size(1)) = model.size(1);
            end
            
            %make a matrix that shows the circular array locations
probe.      plotter     = sparse(probe.indx, probe.indy, true, model.size(1), model.size(2));

            %% Run the FD Model:
            switch solvetype
                case 1 %finite difference 3D
                tic
[press, rfdata] = CenteredDiff3D(model, source, config, probe, pfun);
                timedone = toc;
                case 2 %k-space 
                tic
[press, rfdata] = Kspace(model, source, config, probe, pfun);
                toc
                otherwise
                error("Invalid Solve Type");
            end

disp("Iterations/Sec : " + config.N/timedone)

            %% Save Data to the format that the DAS code expects:
ele         = 512;                                                          % [number of elements]
fs          = 19.2308; %leafdata.fs;                                                 
sos         = config.c/1e3;
x0          = sensarr.x0*1e3;
z0          = sensarr.z0*1e3;
save("ModelRun","ele","fs","rfdata","sos","x0","z0");

            %% Space for testing


            %% Finite Difference function
            function [pressure2, rfdata] = CenteredDiff3D(model, source, config, probe, pfun)
            
            %preallocate an rfdata matrix for the observation timeseries
rfdata      = zeros(config.N,length(probe.indx));
            
            %enforce model stability condition
            if any(model.cfl >= 1/sqrt(3)) 
error("c*dt/dn = " + cfl((cfl >= 1)) + ...
            ", it must be less than or equal to 1/sqrt(3) for numerical stability.")
            end

            if ~isa(pfun, 'function_handle') %input groom
error("pfun must be a valid function handle")
            end

            %get the probe indices
probes      = sub2ind(model.size,probe.indx,probe.indy,probe.indz);

            %iterate over the time-steps
            for ii = 1:config.N
            disp("Running N = "+ii)

            %pause points
%             if mod(ii,100) == 0;
%             pause;
%             end

            % calculate the waves at x,y,z + delta(x,y,z) etc with periodic
            %boundaries
            for dim = 1:3
forward.    (char(119+dim))     = circshift(source.p1,1,dim);               %[Pa] forward-moving wave
backward.   (char(119+dim))     = circshift(source.p1,-1,dim);              %[Pa] backward-moving wave
            end

            if model.PML %if we said to include an absorbing layer for the x and y directions
dpth        = model.Pml.Depth;
            % apply the window to the x and y components of the
            % forward/backward waves:
forward.    x([end:-1:end-dpth+1 1:dpth],:,:) ...
            = forward.x([end:-1:end-dpth+1 1:dpth],:,:).*repmat(model.Pml.halfx,2,1);
backward.   x([end:-1:end-dpth+1 1:dpth],:,:) ...
            = backward.x([end:-1:end-dpth+1 1:dpth],:,:).*repmat(model.Pml.halfx,2,1);
forward.    y(:,[end:-1:end-dpth+1 1:dpth],:) ...
            = forward.y(:,[end:-1:end-dpth+1 1:dpth],:).*repmat(model.Pml.halfy,1,2);
backward.   y(:,[end:-1:end-dpth+1 1:dpth],:) ...
            = backward.y(:,[end:-1:end-dpth+1 1:dpth],:).*repmat(model.Pml.halfy,1,2);
            end
            
            % include the z boundaries
            if model.zboundaries
forward.    z(:,:,1)            = source.p1(:,:,1);                         % hard boundary
backward.   z(:,:,end)          = zeros(size(forward.z,1:2));               % pressure release boundary
            end
            
            % 3D Finite Difference to propagate the wave
pressure2   = 2*source.p1 - source.p0 ...
            + model.cfl(1).^2*(forward.x -2*source.p1 + backward.x)...
            + model.cfl(2).^2*(forward.y -2*source.p1 + backward.y)...
            + model.cfl(3).^2*(forward.z -2*source.p1 + backward.z);

pressure2   = (1-1/1000)*pressure2;                                         % apply a small, artifical attenuation
source.     p0       = source.p1;                                           % increment the pressure field

            if ii <= config.pulse_len                                       % if we're in the pulse length window
pressure2   (source.ind)    = pfun(config.t(ii+1));                         % [Pa] apply the source
            end

source.     p1       = pressure2;

            % Extract pressure data at the observation pointts
rfdata(ii,:)=pressure2(probes);

            %animate the plot in real-time!
            if config.animate
                %round(model.size(3)/2)
imagesc(real(squeeze(source.p1(:,157,:)))'); colormap('pucolors.inferno'); daspect([1 1 1]); clim([-.01 0.01]);
            subtitle("Step " + (ii+1)); title("Hard/Pressure Release Bounds, y = L_y/2"); xlabel('X index'); ylabel('Y index');
            drawnow;
            end
            
            end

            end

            %% K-Space Function
            function [f_s2, rfdata] = Kspace(model, source, config, probe, pfun) 
            %3D K-Space implementation

            if any(model.cfl > 2/pi*config.c/max(source.cxyz,[],"all"))
error("CFL too high!")
            end

            if ~isa(pfun, 'function_handle') %input groom
error("pfun must be a valid function handle")
            end

            % preallocate for rfdata and get probe indices
rfdata      = zeros(config.N,length(probe.indx));
probes      = sub2ind(model.size,probe.indx,probe.indy,probe.indz);

            % pre-iteration setup
rho_root    = sqrt(source.rhoxyz);                                          % convenient variable for density
ccontrast   = (1 - source.cxyz.^2./config.c.^2);                            % speed of sound contrast
rhocontrast = config.c^2.*rho_root.*source.laplacianrho;                    % density contrast
f_s1        = source.p0./rho_root;                                          % normalized scattered field 0
f_s2        = source.p1./rho_root;                                          % normalized scattered field 1
w0          = 0; %ccontrast.*f_s1;                                              % little w
W0          = fftshift(fftn(w0));                                           % big W

            %iterate over time steps:
            for ii = 1:config.N
            disp("Running N = "+ii);

            %increment variables
v           = ccontrast.*f_s2;                                              % little v
q           = rhocontrast.*(f_s2 + w0 - v);                                 % little q
V           = fftn(v);                                                      % big V
Q           = config.c.^2*fftn(q);                                          % big Q
W1          = fftn(f_s2 + v);                                               %calculate W at current timestep

            % Propagate!
W2          = 2*W1 - W0 + 4*sin(config.omg*config.dt/2).^2.* ...            % propagate W with nonstandard second order FD
                (V - W1 - Q./(config.omg.^2));

            % increment again
W0          = W1;                                                           % increment W0
w0          = ifftn(ifftshift(W2));                                         % increment little w0
f_s2        = w0 - v;                                                       % calculate the scattered field

            %extract pressure at observation points (need to denormalize)
rfdata(ii,:)= f_s2(probes).*rho_root(probes);                       

            %apply the new source conditions
            if ii <= config.pulse_len
f_s2         (source.ind)    = pfun(config.t(ii+1))./rho_root(source.ind);
            end

            %animate!
            if config.animate
imagesc(probe.plotter + real(squeeze(f_s2(:,:,floor(model.size(3)/2)))')); colormap('pucolors.inferno'); daspect([1 1 1]); clim([-1 1]); colorbar; drawnow
            end

            end
            
            end

            %% x0 and z0 if you don't have the files

            function [x0, z0] = sensordata()
x0 = [
           -0.0384
           -0.0385
           -0.0384
           -0.0389
           -0.0383
           -0.0383
           -0.0388
           -0.0387
           -0.0387
           -0.0386
           -0.0386
           -0.0385
           -0.0384
           -0.0383
           -0.0380
           -0.0379
           -0.0378
           -0.0379
           -0.0378
           -0.0377
           -0.0375
           -0.0374
           -0.0370
           -0.0372
           -0.0370
           -0.0368
           -0.0367
           -0.0365
           -0.0363
           -0.0361
           -0.0360
           -0.0358
           -0.0356
           -0.0354
           -0.0352
           -0.0350
           -0.0348
           -0.0345
           -0.0344
           -0.0341
           -0.0339
           -0.0337
           -0.0334
           -0.0332
           -0.0329
           -0.0326
           -0.0324
           -0.0321
           -0.0318
           -0.0313
           -0.0313
           -0.0310
           -0.0308
           -0.0305
           -0.0302
           -0.0298
           -0.0295
           -0.0292
           -0.0290
           -0.0286
           -0.0283
           -0.0280
           -0.0277
           -0.0273
           -0.0270
           -0.0266
           -0.0263
           -0.0259
           -0.0256
           -0.0252
           -0.0248
           -0.0245
           -0.0242
           -0.0238
           -0.0234
           -0.0230
           -0.0226
           -0.0223
           -0.0219
           -0.0215
           -0.0211
           -0.0207
           -0.0203
           -0.0199
           -0.0195
           -0.0191
           -0.0187
           -0.0182
           -0.0178
           -0.0174
           -0.0170
           -0.0166
           -0.0161
           -0.0157
           -0.0153
           -0.0148
           -0.0144
           -0.0140
           -0.0135
           -0.0131
           -0.0127
           -0.0122
           -0.0118
           -0.0113
           -0.0109
           -0.0104
           -0.0100
           -0.0095
           -0.0091
           -0.0086
           -0.0081
           -0.0077
           -0.0072
           -0.0068
           -0.0063
           -0.0058
           -0.0054
           -0.0049
           -0.0044
           -0.0040
           -0.0035
           -0.0030
           -0.0026
           -0.0021
           -0.0016
           -0.0012
           -0.0007
           -0.0002
            0.0002
            0.0007
            0.0012
            0.0016
            0.0021
            0.0026
            0.0030
            0.0035
            0.0040
            0.0044
            0.0049
            0.0054
            0.0058
            0.0063
            0.0068
            0.0072
            0.0077
            0.0081
            0.0086
            0.0091
            0.0095
            0.0100
            0.0104
            0.0109
            0.0113
            0.0118
            0.0122
            0.0127
            0.0131
            0.0135
            0.0140
            0.0144
            0.0148
            0.0153
            0.0157
            0.0161
            0.0166
            0.0170
            0.0174
            0.0178
            0.0182
            0.0186
            0.0191
            0.0195
            0.0199
            0.0203
            0.0207
            0.0210
            0.0214
            0.0218
            0.0222
            0.0226
            0.0230
            0.0234
            0.0238
            0.0241
            0.0245
            0.0249
            0.0252
            0.0256
            0.0259
            0.0262
            0.0266
            0.0269
            0.0273
            0.0276
            0.0280
            0.0283
            0.0286
            0.0289
            0.0292
            0.0295
            0.0299
            0.0302
            0.0305
            0.0307
            0.0311
            0.0313
            0.0316
            0.0319
            0.0322
            0.0324
            0.0326
            0.0329
            0.0332
            0.0334
            0.0337
            0.0339
            0.0341
            0.0343
            0.0346
            0.0348
            0.0350
            0.0352
            0.0355
            0.0356
            0.0358
            0.0360
            0.0362
            0.0364
            0.0366
            0.0367
            0.0368
            0.0370
            0.0372
            0.0373
            0.0375
            0.0375
            0.0377
            0.0378
            0.0379
            0.0381
            0.0381
            0.0382
            0.0383
            0.0384
            0.0385
            0.0386
            0.0386
            0.0387
            0.0388
            0.0389
            0.0388
            0.0389
            0.0390
            0.0390
            0.0390
            0.0390
            0.0391
            0.0391
            0.0390
            0.0390
            0.0389
            0.0389
            0.0389
            0.0388
            0.0388
            0.0387
            0.0386
            0.0385
            0.0385
            0.0384
            0.0383
            0.0382
            0.0381
            0.0380
            0.0379
            0.0377
            0.0376
            0.0375
            0.0374
            0.0373
            0.0370
            0.0369
            0.0368
            0.0366
            0.0364
            0.0363
            0.0361
            0.0359
            0.0357
            0.0355
            0.0353
            0.0351
            0.0349
            0.0346
            0.0344
            0.0342
            0.0340
            0.0338
            0.0335
            0.0333
            0.0330
            0.0328
            0.0325
            0.0323
            0.0320
            0.0317
            0.0315
            0.0312
            0.0308
            0.0306
            0.0303
            0.0300
            0.0297
            0.0294
            0.0291
            0.0288
            0.0285
            0.0281
            0.0278
            0.0274
            0.0271
            0.0268
            0.0264
            0.0261
            0.0257
            0.0254
            0.0250
            0.0246
            0.0243
            0.0239
            0.0235
            0.0232
            0.0228
            0.0224
            0.0220
            0.0216
            0.0212
            0.0209
            0.0205
            0.0201
            0.0197
            0.0192
            0.0188
            0.0184
            0.0180
            0.0176
            0.0172
            0.0168
            0.0163
            0.0159
            0.0154
            0.0150
            0.0146
            0.0141
            0.0137
            0.0133
            0.0128
            0.0124
            0.0119
            0.0115
            0.0110
            0.0106
            0.0101
            0.0097
            0.0092
            0.0088
            0.0083
            0.0079
            0.0074
            0.0069
            0.0065
            0.0060
            0.0056
            0.0051
            0.0046
            0.0042
            0.0037
            0.0032
            0.0028
            0.0023
            0.0018
            0.0014
            0.0009
            0.0004
           -0.0001
           -0.0005
           -0.0010
           -0.0015
           -0.0019
           -0.0024
           -0.0029
           -0.0033
           -0.0038
           -0.0043
           -0.0047
           -0.0052
           -0.0057
           -0.0061
           -0.0066
           -0.0070
           -0.0075
           -0.0080
           -0.0084
           -0.0089
           -0.0093
           -0.0098
           -0.0102
           -0.0107
           -0.0111
           -0.0116
           -0.0120
           -0.0125
           -0.0129
           -0.0134
           -0.0138
           -0.0142
           -0.0147
           -0.0151
           -0.0155
           -0.0159
           -0.0164
           -0.0168
           -0.0172
           -0.0176
           -0.0181
           -0.0185
           -0.0189
           -0.0193
           -0.0197
           -0.0201
           -0.0205
           -0.0209
           -0.0213
           -0.0217
           -0.0220
           -0.0224
           -0.0228
           -0.0232
           -0.0236
           -0.0239
           -0.0243
           -0.0247
           -0.0250
           -0.0254
           -0.0257
           -0.0261
           -0.0264
           -0.0268
           -0.0271
           -0.0274
           -0.0277
           -0.0281
           -0.0284
           -0.0287
           -0.0290
           -0.0293
           -0.0296
           -0.0299
           -0.0302
           -0.0306
           -0.0308
           -0.0311
           -0.0314
           -0.0316
           -0.0320
           -0.0322
           -0.0325
           -0.0327
           -0.0329
           -0.0332
           -0.0334
           -0.0337
           -0.0339
           -0.0341
           -0.0343
           -0.0345
           -0.0347
           -0.0350
           -0.0352
           -0.0354
           -0.0355
           -0.0357
           -0.0359
           -0.0361
           -0.0363
           -0.0364
           -0.0366
           -0.0367
           -0.0369
           -0.0371
           -0.0372
           -0.0373
           -0.0374
           -0.0375
           -0.0376
           -0.0378
           -0.0379
           -0.0380
           -0.0380
           -0.0381
           -0.0382
           -0.0383
           -0.0383
           -0.0384
           -0.0385
           -0.0385
           -0.0386
           -0.0386
           -0.0386
           -0.0386
           -0.0387
           -0.0387  ].';

z0 = [
            0.0015
            0.0019
            0.0024
            0.0029
            0.0033
            0.0038
            0.0043
            0.0048
            0.0052
            0.0057
            0.0062
            0.0066
            0.0071
            0.0075
            0.0080
            0.0084
            0.0089
            0.0094
            0.0098
            0.0103
            0.0107
            0.0112
            0.0116
            0.0121
            0.0125
            0.0130
            0.0134
            0.0138
            0.0143
            0.0147
            0.0152
            0.0156
            0.0160
            0.0164
            0.0169
            0.0173
            0.0177
            0.0181
            0.0186
            0.0190
            0.0194
            0.0198
            0.0202
            0.0206
            0.0210
            0.0213
            0.0217
            0.0221
            0.0225
            0.0227
            0.0233
            0.0237
            0.0240
            0.0244
            0.0247
            0.0251
            0.0255
            0.0258
            0.0262
            0.0265
            0.0269
            0.0272
            0.0275
            0.0279
            0.0282
            0.0285
            0.0288
            0.0291
            0.0294
            0.0297
            0.0300
            0.0303
            0.0307
            0.0310
            0.0312
            0.0315
            0.0318
            0.0321
            0.0323
            0.0326
            0.0328
            0.0331
            0.0333
            0.0336
            0.0338
            0.0340
            0.0342
            0.0345
            0.0347
            0.0349
            0.0351
            0.0353
            0.0355
            0.0357
            0.0359
            0.0360
            0.0362
            0.0364
            0.0366
            0.0367
            0.0369
            0.0370
            0.0372
            0.0373
            0.0374
            0.0376
            0.0377
            0.0378
            0.0379
            0.0380
            0.0381
            0.0382
            0.0383
            0.0384
            0.0385
            0.0386
            0.0386
            0.0387
            0.0387
            0.0388
            0.0388
            0.0388
            0.0389
            0.0389
            0.0389
            0.0389
            0.0389
            0.0390
            0.0390
            0.0390
            0.0390
            0.0389
            0.0389
            0.0389
            0.0388
            0.0388
            0.0387
            0.0387
            0.0386
            0.0386
            0.0385
            0.0384
            0.0384
            0.0383
            0.0382
            0.0381
            0.0380
            0.0379
            0.0378
            0.0376
            0.0375
            0.0374
            0.0373
            0.0371
            0.0370
            0.0369
            0.0367
            0.0366
            0.0363
            0.0362
            0.0360
            0.0358
            0.0357
            0.0355
            0.0353
            0.0351
            0.0348
            0.0346
            0.0344
            0.0342
            0.0340
            0.0338
            0.0335
            0.0333
            0.0331
            0.0328
            0.0325
            0.0323
            0.0320
            0.0318
            0.0315
            0.0312
            0.0309
            0.0306
            0.0303
            0.0300
            0.0297
            0.0294
            0.0291
            0.0288
            0.0285
            0.0282
            0.0278
            0.0275
            0.0272
            0.0269
            0.0265
            0.0261
            0.0258
            0.0255
            0.0251
            0.0248
            0.0244
            0.0240
            0.0237
            0.0233
            0.0229
            0.0225
            0.0221
            0.0217
            0.0213
            0.0210
            0.0206
            0.0202
            0.0198
            0.0194
            0.0189
            0.0185
            0.0181
            0.0177
            0.0173
            0.0169
            0.0165
            0.0160
            0.0156
            0.0152
            0.0147
            0.0143
            0.0139
            0.0134
            0.0130
            0.0125
            0.0121
            0.0117
            0.0112
            0.0107
            0.0103
            0.0098
            0.0094
            0.0089
            0.0085
            0.0080
            0.0075
            0.0071
            0.0066
            0.0062
            0.0057
            0.0052
            0.0048
            0.0043
            0.0038
            0.0034
            0.0029
            0.0024
            0.0020
            0.0015
           -0.0015
           -0.0020
           -0.0025
           -0.0029
           -0.0034
           -0.0039
           -0.0043
           -0.0048
           -0.0053
           -0.0057
           -0.0062
           -0.0067
           -0.0071
           -0.0076
           -0.0080
           -0.0085
           -0.0090
           -0.0094
           -0.0099
           -0.0103
           -0.0108
           -0.0112
           -0.0117
           -0.0121
           -0.0125
           -0.0130
           -0.0134
           -0.0139
           -0.0143
           -0.0148
           -0.0152
           -0.0156
           -0.0160
           -0.0165
           -0.0169
           -0.0173
           -0.0177
           -0.0181
           -0.0185
           -0.0190
           -0.0194
           -0.0198
           -0.0202
           -0.0206
           -0.0210
           -0.0214
           -0.0217
           -0.0221
           -0.0225
           -0.0229
           -0.0233
           -0.0237
           -0.0240
           -0.0244
           -0.0247
           -0.0251
           -0.0255
           -0.0258
           -0.0262
           -0.0265
           -0.0269
           -0.0272
           -0.0275
           -0.0279
           -0.0282
           -0.0285
           -0.0288
           -0.0292
           -0.0295
           -0.0298
           -0.0301
           -0.0303
           -0.0306
           -0.0309
           -0.0312
           -0.0315
           -0.0318
           -0.0320
           -0.0323
           -0.0326
           -0.0328
           -0.0331
           -0.0333
           -0.0336
           -0.0338
           -0.0341
           -0.0343
           -0.0345
           -0.0347
           -0.0350
           -0.0352
           -0.0354
           -0.0355
           -0.0357
           -0.0359
           -0.0361
           -0.0362
           -0.0364
           -0.0366
           -0.0367
           -0.0369
           -0.0371
           -0.0372
           -0.0373
           -0.0375
           -0.0376
           -0.0377
           -0.0378
           -0.0380
           -0.0381
           -0.0382
           -0.0383
           -0.0384
           -0.0384
           -0.0385
           -0.0387
           -0.0387
           -0.0388
           -0.0388
           -0.0389
           -0.0389
           -0.0390
           -0.0390
           -0.0390
           -0.0391
           -0.0391
           -0.0391
           -0.0391
           -0.0390
           -0.0391
           -0.0390
           -0.0390
           -0.0390
           -0.0389
           -0.0389
           -0.0388
           -0.0388
           -0.0388
           -0.0387
           -0.0387
           -0.0386
           -0.0385
           -0.0385
           -0.0384
           -0.0383
           -0.0382
           -0.0381
           -0.0380
           -0.0378
           -0.0377
           -0.0376
           -0.0375
           -0.0373
           -0.0372
           -0.0370
           -0.0369
           -0.0368
           -0.0366
           -0.0364
           -0.0363
           -0.0361
           -0.0359
           -0.0357
           -0.0355
           -0.0353
           -0.0352
           -0.0349
           -0.0347
           -0.0345
           -0.0343
           -0.0340
           -0.0338
           -0.0336
           -0.0334
           -0.0331
           -0.0329
           -0.0326
           -0.0323
           -0.0321
           -0.0318
           -0.0315
           -0.0312
           -0.0310
           -0.0307
           -0.0304
           -0.0301
           -0.0298
           -0.0295
           -0.0291
           -0.0289
           -0.0285
           -0.0282
           -0.0279
           -0.0275
           -0.0272
           -0.0269
           -0.0266
           -0.0262
           -0.0259
           -0.0255
           -0.0251
           -0.0248
           -0.0244
           -0.0241
           -0.0237
           -0.0233
           -0.0229
           -0.0225
           -0.0222
           -0.0218
           -0.0214
           -0.0210
           -0.0206
           -0.0202
           -0.0198
           -0.0194
           -0.0190
           -0.0186
           -0.0181
           -0.0177
           -0.0173
           -0.0169
           -0.0165
           -0.0160
           -0.0156
           -0.0152
           -0.0147
           -0.0143
           -0.0139
           -0.0134
           -0.0130
           -0.0126
           -0.0121
           -0.0117
           -0.0112
           -0.0108
           -0.0103
           -0.0099
           -0.0094
           -0.0090
           -0.0085
           -0.0080
           -0.0076
           -0.0071
           -0.0067
           -0.0062
           -0.0057
           -0.0053
           -0.0048
           -0.0043
           -0.0039
           -0.0034
           -0.0029
           -0.0025
           -0.0020
           -0.0015  ].';
            end
