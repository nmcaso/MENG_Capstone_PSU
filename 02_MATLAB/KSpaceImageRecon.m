            %% Setup
            clear; clc;
            for folder = strtrim(string(ls('*0*'))).'
                addpath(folder);
            end

            %% Initiate Classes
leafdata    = Dataset("crossinghair",0);
sensarr     = SensorArray("crossinghair");
mins        = -0.02; intval = 0.0001; maxs = -mins-intval;
area1       = ImageArea(maxs, mins, intval, maxs, mins, intval);
ind_mat     = IndexMatrix(sensarr,area1,leafdata, "index", 0); disp("Time to make Index Matrix: " + ind_mat.Times.total)

            %% Setup
c           = cast(leafdata.c,'single');
dt          = cast(leafdata.dt,'single');
s1          = size(leafdata.rfdata,1)*2;
s2          = size(leafdata.rfdata,2);

[theta, r]  = cart2pol(sensarr.x0, sensarr.z0);
[theta, s]  = sort(theta);
r           = mean(r);

k_time      = c*dt*(-s1/2:s1/2-1);
k_time      = repmat(k_time(:),[1 s2]);
k_theta     = repmat(theta,[s1 1]);

k           = hypot( k_theta, k_time);

omega       = leafdata.c .* k_time;
omega_int   = leafdata.c .* k;

scale       = leafdata.c.^2 .* sqrt((omega./leafdata.c).^2 - k_theta.^2)./(2.*omega);
scale       (omega == 0 & k_theta) = leafdata.c./2;

indicesvec  = 1:(s1)/2;
thetas      = repmat(theta,[(s1)/2 1]);


            %% Vectorized algorithm
tic
rfpad       = [flip(leafdata.rfdata,1); leafdata.rfdata];
RFDATA      = fftshift(nufftn(ifftshift(rfpad),omega(:),omega_int(:)));

KSPCIMG     = real(fftshift(ifftn(ifftshift(RFINT))))*4./leafdata.c;
KSPCIMG     = KSPCIMG(indicesvec,:);
toc

