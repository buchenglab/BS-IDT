clc
clearvars
close all

addpath('functions')
addpath('Utilities')

%% Handle Declarations
FT = @(x) fftshift(fft2(ifftshift(x)));
iFT = @(x) fftshift(ifft2(ifftshift(x)));

%% Parameter

lambda=0.447; % Wavelength
k=2*pi/lambda; % Wave number

NA=0.655;
Mag=40;
Pixelsize=6.5/Mag;
n_Medium=1.34;
dz = lambda./(NA.^2);

Length_MN=32;% the number of LEDs

Calib= 0; % if Calibrate LED postion
gpu = 0;  % if using gpu
Bright_Radius=10;% the Radius of LED ring
fDL=Bright_Radius./tan(asin(NA));% the distance of LED and object 

% Sorted_Pos = [360:-22.5:1];
% Sorted_Pos = repelem(Sorted_Pos, 2);
% Sorted_Pos = [cosd(Sorted_Pos'), sind(Sorted_Pos')];
% Sorted_Pos = Sorted_Pos';
%% Step Load measured intensity data and initial spectrum postion

genfol = 'K:\BU\Data\LIDT\BC_Mitotic\';  % 'D:\BU\Data\LIDT\Rich_Lipid_img\';  %  
prefix = '210424_Mitotic_1';  % '210609_BCancerImg';
suffix = 'acq-crop-precal';
avg = 60;
wvlngth = {'noIR'};
illum_idx = [1:32];

% Load data
tmpnm = [genfol, prefix,'_', wvlngth{1}, '_av', num2str(avg),'_',suffix];
tmpfol = dir(tmpnm);
dat_fol = [tmpnm, '\', tmpfol(3).name, '\'];
load([dat_fol 'IDT_Data_pls600.mat'], 'img', 'iNA');

I_Raw = img(:,:,illum_idx);

% Normalize image without background subtraction
I_Raw = BkgndNorm(I_Raw, false) + 1; 

% Crop image down
cent = [996, 1302];  %[855, 1015]; %[801, 801]; %[1000, 1312]; %[1164, 1356];%%[1028, 1293]; %Center:  %Size: 2160 x 2560
sz = [1024, 1024];
I_Raw = I_Raw(cent(1) - sz(1)/2:cent(1) + sz(1)/2 - 1,...
          cent(2) - sz(2)/2:cent(2) + sz(2)/2 - 1,:); 

% Convert iNA to spatial frequencies, correct 
iNA = -iNA/lambda;

Cablib_Nx=256;
Cablib_Ny=256;

Cablib_pointX=513;  
Cablib_pointY=513;
% Sorted_Pos = Sorted_Pos';
eval Step0_IDT_Init_v2

%% Step Implete the calibation of LED postion

 if Calib==1
    eval Step1_IDT_Calib
 end
 
%% Implete the IDT
Alpha=[1];

nAxes = 2;  % Select number of asymmetry axes to use
dat_fol = [dat_fol 'Reprocessed_Results\'];
mkdir(dat_fol);
for nA = 1:length(Alpha)
    disp(['Processing ' num2str(Alpha(nA)) ' Regularization...']);
    % Reassign angles, Illum quantity to prevent shrinkage in getDPC script
    Length_MN=32;
    Ini_NAx = iNA(:,1);
    Ini_NAy = iNA(:,2);
    
    % Process hot and cold measurements, take difference
    eval Step2_GetDPC_v2

    %% Save Results

    hotfile = Tiff([dat_fol 'DPC_inv_hot_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.tiff'],'w');
    writeTiff_32bGray(hotfile, V_hot);
    coldfile = Tiff([dat_fol 'DPC_inv_cold_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.tiff'],'w');
    writeTiff_32bGray(coldfile, V_cold);
    difffile = Tiff([dat_fol 'DPC_inv_diff_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.tiff'],'w');
    writeTiff_32bGray(difffile, V_diff);
    hotfile = Tiff([dat_fol 'DPC_inv_hot_RI_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.tiff'],'w');
    writeTiff_32bGray(hotfile, RI_hot);
    coldfile = Tiff([dat_fol 'DPC_inv_cold_RI_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.tiff'],'w');
    writeTiff_32bGray(coldfile, RI_cold);
    difffile = Tiff([dat_fol 'DPC_inv_diff_RI_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.tiff'],'w');
    writeTiff_32bGray(difffile, diffRI);

    save([dat_fol 'Recon_DPC_inv_Tau' num2str(Alpha(nA)) '_cal' num2str(Calib) '.mat'],'V_hot','V_cold','V_diff','RI_hot','RI_cold','diffRI','Calib');
end


