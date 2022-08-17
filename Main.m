%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%           Bond-Selective Intensity Diffraction Tomography Script
%
%
%   Purpose: This is an example reconstruction script for processing data
%            obtained from the BS-IDT system. This code accepts "hot" and
%            "cold" raw intensity images at multiple illumination angles
%            where the IR pump beam is on and off, respectively, and
%            reconstructs a "hot" and "cold" refractive index (RI) volume
%            of the sample. The difference of these volumes is then
%            taken to recover the chemical signal from the object of
%            interest. This code has been adapted primarily from the 
%            annular IDT script originally published in the work 
%            "High speed in vitro intensity diffraction tomography" by Li
%            and Matlock et al. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

%% Add functions and utilities necessary for BS-IDT to the path
addpath('functions')
addpath('Utilities')
addpath('data')

%% Handle Declarations
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

%% Variable Intialization

% Microscope parameters
lambda=0.447; % Set imaging wavelength (um)
NA=0.65;  % Microscope objective NA
Mag=40;   % Microscope magnification
Pixelsize=6.5/Mag;  % Set pixel size at object plane (um)
Bright_Radius=10;% LED ring radius

% Imaging medium, image parameters
n_Medium=1.34;  % Assumed imaging medium refractive index
sz = [1024, 1024]; % set reconstructed image size

k=2*pi/lambda; % Set image wavenumber
fDL=Bright_Radius./tan(asin(NA));% the distance of LED and object 

% File name for processing
fnm = 'Cells_1657cm-1';

% Regularization values for IDT reconstruction
Tau = [1e2, 1e2]; % Real and imaginary regularization parameters

% set toggles for calibrating illumination angle or using gpu
Calib= 0; % 1 if calibrating the LED position, 0 if not
gpu = 0;  % 1 if using gpu, 0 if not

%% Step 1: Load measured intensity data and initial spectrum postion

  % Set file name for data to process

load([cd '\data\' fnm '.mat'], 'img', 'iNA', 'Sorted_Pos');
cent = size(img)/2+1;  % Calculate image center
Length_MN = size(iNA, 1);  % Obtain total number of illuminations

% Crop raw images to desired reconstruction size
I_Raw = img(cent(1) - sz(1)/2:cent(1) + sz(1)/2 - 1,...
          cent(2) - sz(2)/2:cent(2) + sz(2)/2 - 1,:); 

% Remove image background
I_Raw = BkgndNorm(I_Raw, false); 
clear img cent
 
% Change illumination angles to work with function
iNA = -iNA/lambda;
if(Calib == 0)
    iNA = [iNA(:,1), -iNA(:,2)];
end

% Set window size for illumination angle calibration
Calib_Nx=sz(1);  
Calib_Ny=sz(2);
Calib_pointX=sz(1)/2+1;  
Calib_pointY=sz(2)/2+1;

eval Step1_IDT_Init

%% Step 2: Implement LED calibration

 if Calib==1
    eval Step2_IDT_Calib
 end
 
%% Step 3: Perform IDT reconstructions
dz=(lambda/(NA^2))/2;
Depth_Set=[-10:10] * dz/n_Medium;

% Process hot and cold measurements, take difference
eval Step3_IDT_Poss

% Take difference of hot, cold measurements
RI_diff = real(RI_hot) - real(RI_cold);

%% Step 4: Save Results
save(['Reconstruction.mat'],'RI_hot','RI_cold','RI_diff','Calib','-v7.3');
   

