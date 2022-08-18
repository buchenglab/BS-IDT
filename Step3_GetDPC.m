disp('Step 3: Applying DPC reconstruction...')

% Select Calibrated illumination NAs
 if Calib==1
    Ini_NAy=(freqXY3(:,1)-(Calib_Nx/2+1))/Calib_Nx/Pixelsize;
    Ini_NAx=(freqXY3(:,2)-(Calib_Ny/2+1))/Calib_Ny/Pixelsize;
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
   
    tmp = conj(F(source_m .* Aperture_fun)).* F(Aperture_fun);
    DC(i) = sum(sum(source_m.*abs(Aperture_fun).^2));
    Hi(:,:,i) = 2 * Ft(real(tmp))/DC(i);
    Hr(:,:,i) = 1i*2 * Ft(1i * imag(tmp))/DC(i);

end


%% Sum images, TFs for DPC image generation
[I_Hot, idx, ~, ~] = applyDPC_Gen((I_Raw(:,:,1:2:32)), [Ini_NAx, Ini_NAy] .* lambda, [0.6, 0.7], pi/nAxes);
[I_Cold, idx, ~, ~] = applyDPC_Gen((I_Raw(:,:,2:2:32)), [Ini_NAx, Ini_NAy] .* lambda, [0.6, 0.7], pi/nAxes);

TF = zeros(size(I_Raw,1), size(I_Raw,2), 4);
for k = 1:size(idx,2)
    tmp_idx = idx(:,k);
    tmp_idx(tmp_idx == 0) = [];
    TF(:,:,k) = sum(Hr(:,:,tmp_idx)./sum(DC(tmp_idx)), 3);
end

I_DPC_Hot = F(I_Hot);
I_DPC_Cold = F(I_Cold);

cnt = 1;
for k = 1:2:size(idx,2)
    DPC_TF(:,:,cnt) = TF(:,:,k) - TF(:,:,k+1);  % (TF(:,:,1) - TF(:,:,2));
    cnt = cnt + 1;
end

Phi_hot = sum(conj(DPC_TF) .* I_DPC_Hot,3) ./ (sum(abs(DPC_TF).^2,3) + Tau);
Phi_cold = sum(conj(DPC_TF) .* I_DPC_Cold,3) ./ (sum(abs(DPC_TF).^2,3) + Tau);

Phi_hot = real(Ft(Phi_hot));
Phi_cold = real(Ft(Phi_cold));
Phi_diff = Phi_hot - Phi_cold;
