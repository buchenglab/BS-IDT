
%% Load measured intensity data and define region for calibration
disp('Step 1: Generating sampling grids and illumination angles...');
imSz=size(I_Raw);
Nx=imSz(1);
Ny=imSz(2);
Nz=Nx;

I_Calib=I_Raw(Calib_pointX-Calib_Nx/2:Calib_pointX+Calib_Nx/2-1,Calib_pointY-Calib_Ny/2:Calib_pointY+Calib_Ny/2-1,:);
  
%% Calculate Fourier space grids, points on Ewald sphere surface

% Frequency coordinate
Max_frequency=NA/lambda;
delta_x = 1/(Pixelsize*Nx);      % frequency sampling X.
delta_y = 1/(Pixelsize*Ny);      % frequency sampling Y.
delta_z = 1/(Pixelsize*Nz);      % frequency sampling Z.

fx = (-fix(Nx/2):1:fix((Nx-1)/2))*delta_x; % frequency coordinate X.
fy = (-fix(Ny/2):1:fix((Ny-1)/2))*delta_y; % frequency coordinate Y.
[fx2D, fy2D] = meshgrid(fx,fy);

% Generating coordinates on the surface of Ewald Sphere
fz2D = real(sqrt((n_Medium./lambda).^2-fy2D.^2-fx2D.^2));
[Theta,R] = cart2pol(fx2D,fy2D);
Aperture = ~(R>Max_frequency);
Aperture_fun = double(Aperture);

%% Calculate illumination initial positions

%initial guess of led postion and frequecy coord
Ini_PixelShiftx=zeros(1,Length_MN);
Ini_PixelShifty=zeros(1,Length_MN);

Ini_NAx = iNA(:,1);
Ini_NAy = iNA(:,2);
Ini_NAz = real(sqrt((n_Medium/lambda)^2 - Ini_NAx.^2 - Ini_NAy.^2));


