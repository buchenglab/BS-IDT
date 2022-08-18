
%%%%% IDT poss after LED Postion calibtion %%%%%

% Select Calibrated illumination NAs
 if Calib==1
    Ini_NAx=(freqXY3(:,1)'-(Cablib_Nx/2+1))/Cablib_Nx/Pixelsize;
    Ini_NAy=(freqXY3(:,2)'-(Cablib_Ny/2+1))/Cablib_Ny/Pixelsize;
 end

% Select only unique illuminations, remove repetitions
Ini_NAx = Ini_NAx(1:2:32);
Ini_NAy = Ini_NAy(1:2:32);
Length_MN = Length_MN/2;

%%
% Convert illumination NA to nearest pixel position
Ini_PixelShiftx=zeros(1,Length_MN);
Ini_PixelShifty=zeros(1,Length_MN);

for i=1:Length_MN
    pic_pos=i;
    
    Ini_PixelShiftx(pic_pos)=round(Ini_NAx(pic_pos)*Pixelsize*Nx);
    Ini_PixelShifty(pic_pos)=round(Ini_NAy(pic_pos)*Pixelsize*Ny);
end

%% Initialize Transfer Functions
du = 1/(Nx * Pixelsize);

% Set source radius, generate mask
s_rad = du/2; %should be 5mm, trying first with single pixel source

uu=fx2D;
vv=fy2D;
k_Wavenum=k;
k_Medium=k_Wavenum*n_Medium;

source = double(real(sqrt((uu./s_rad).^2 + (vv./s_rad).^2)) <= 1);

Hr = zeros(Ny, Nx, Length_MN);
Hi = zeros(Ny, Nx, Length_MN);

for i=1:Length_MN
    pic_pos=i;

    v_s=Ini_NAx(pic_pos);
    u_s=Ini_NAy(pic_pos);
    
    source_p = circshift(source,[Ini_PixelShifty(pic_pos),Ini_PixelShiftx(pic_pos)]);
    source_m = circshift(source,-[Ini_PixelShifty(pic_pos),Ini_PixelShiftx(pic_pos)]);
   
    tmp = conj(FT(source_m .* Aperture_fun)).* FT(Aperture_fun);
    DC(i) = sum(sum(source_m.*abs(Aperture_fun).^2));
    Hi(:,:,i) = 2 * iFT(real(tmp))/DC(i);
    Hr(:,:,i) = 1i*2 * iFT(1i * imag(tmp))/DC(i);
%     Hr = 2*Ft(real(FSPc.*FP));
% Hi = -2*Ft(imag(FSPc.*FP));
% 
% Htot = sqrt(abs(Hr).^2+abs(Hi).^2);
% Htotmax = max(max(Htot));
% Hr = Hr./Htotmax;
% Hi = Hi./Htotmax;
%     disp(['NaN: ' num2str(sum(isnan(reshape(abs(Hr(:,:,i)), [numel(Hr(:,:,i)), 1]))))]);
end


%% Sum images, TFs for DPC image generation
% [I_Hot, idx, ~, ~] = applyDPC_Gen((I_Raw(:,:,1:2:32)), [Ini_NAx, Ini_NAy] .* lambda, [0.6, 0.7], pi/nAxes);
% [I_Cold, idx, ~, ~] = applyDPC_Gen((I_Raw(:,:,2:2:32)), [Ini_NAx, Ini_NAy] .* lambda, [0.6, 0.7], pi/nAxes);

idx = [0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 12, 0, 0, 0, 0;
       1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14, 15, 16;
       0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 11, 12, 13, 14, 15, 16;
       0, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0]';

tmpHot = I_Raw(:,:,1:2:32);
tmpCold = I_Raw(:,:,2:2:32);

tmpTop = sum(tmpHot(:,:,logical(idx(:,1))), 3);
tmpBot = sum(tmpHot(:,:,logical(idx(:,2))), 3);
I_Hot(:,:,1) = (tmpTop - tmpBot)./(tmpTop + tmpBot);
tmpTop = sum(tmpHot(:,:,logical(idx(:,3))), 3);
tmpBot = sum(tmpHot(:,:,logical(idx(:,4))), 3);
I_Hot(:,:,2) = (tmpTop - tmpBot)./(tmpTop + tmpBot);

tmpTop = sum(tmpCold(:,:,logical(idx(:,1))), 3);
tmpBot = sum(tmpCold(:,:,logical(idx(:,2))), 3);
I_Cold(:,:,1) = (tmpTop - tmpBot)./(tmpTop + tmpBot);
tmpTop = sum(tmpCold(:,:,logical(idx(:,3))), 3);
tmpBot = sum(tmpCold(:,:,logical(idx(:,4))), 3);
I_Cold(:,:,2) = (tmpTop - tmpBot)./(tmpTop + tmpBot);


TF = zeros(size(I_Raw,1), size(I_Raw,2), 4);
for k = 1:size(idx,2)
    tmp_idx = idx(:,k);
    tmp_idx(tmp_idx == 0) = [];
    TF(:,:,k) = sum(Hr(:,:,tmp_idx)./sum(DC(tmp_idx)), 3);%
end

% tmp = zeros(1024);
% for k = 1:7
%    tmp = tmp + Hr(:,:,idx(k,3)); 
%    figure(534545)
%    imagesc(imag(tmp));axis image off;
%    pause(0.1);
% end

I_DPC_Hot = FT(I_Hot);
I_DPC_Cold = FT(I_Cold);

cnt = 1;
for k = 1:2:size(idx,2)
    DPC_TF(:,:,cnt) = TF(:,:,k) - TF(:,:,k+1);  % (TF(:,:,1) - TF(:,:,2));
    cnt = cnt + 1;
end

V_hot = sum(conj(DPC_TF) .* I_DPC_Hot,3) ./ (sum(abs(DPC_TF).^2,3) + Alpha(nA));
V_cold = sum(conj(DPC_TF) .* I_DPC_Cold,3) ./ (sum(abs(DPC_TF).^2,3) + Alpha(nA));

V_hot = real(iFT(V_hot));
V_cold = real(iFT(V_cold));
V_diff = V_hot - V_cold;

RI_hot = (V_hot * lambda)./(2 * pi * dz) + n_Medium;
RI_cold = (V_cold * lambda)./(2 * pi * dz) + n_Medium;
diffRI = RI_hot - RI_cold;

imstruct.ImageLength = size(Hr, 1);
imstruct.ImageWidth = size(Hr, 2);
imstruct.Photometric = Tiff.Photometric.LinearRaw;
imstruct.Compression = Tiff.Compression.None;
imstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
imstruct.BitsPerSample = 32;
imstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;

mkdir([dat_fol 'DPC TF_Ax' num2str(nAxes)]);
for q = 1:size(DPC_TF,3)
    freTF = Tiff([dat_fol 'DPC TF_Ax' num2str(nAxes) '\TF_' num2str(q) '.tiff'],'w');
    setTag(freTF, imstruct);
    write(freTF, single(imag(DPC_TF(:,:,q))));
    close(freTF);
end

