
%%%%% IDT poss after LED Postion calibtion %%%%%
disp('Step 3: Generating Transfer Functions and Inverting...');
% apply calibrated
 if Calib==1
    Ini_NAx=(freqXY3(:,1)-(Calib_Nx/2+1))/Calib_Nx/Pixelsize;
    Ini_NAy=(freqXY3(:,2)-(Calib_Ny/2+1))/Calib_Ny/Pixelsize;
 else
    TempNA=Ini_NAx;
    Ini_NAx=-Ini_NAy;
    Ini_NAy=TempNA;
 end

Ini_NAx = Ini_NAx(1:2:size(Ini_NAx,1));
Ini_NAy = Ini_NAy(1:2:size(Ini_NAy,1));
Length_MN = Length_MN/2;

%% Generate Complex or pure phase 3D object, and correspoding intensity bright field iamges
uu=fx2D;
vv=fy2D;
k_Wavenum=k;
k_Medium=k_Wavenum*n_Medium;

disp('Calculating slice-wise transfer functions...');

% Preallocate Transfer function matrices
PTF_4D=single(zeros(Nx,Ny,length(Depth_Set),Length_MN));
ATF_4D=single(zeros(Nx,Ny,length(Depth_Set),Length_MN));

% Preallocate vectors containing pixel shifts based on illuminations
Ini_PixelShiftx=zeros(1,Length_MN);
Ini_PixelShifty=zeros(1,Length_MN);

for nL=1:Length_MN

    % Convert calibrated illumination spatial freq. to pixel values
    Ini_PixelShiftx(nL)=round(Ini_NAx(nL)*Pixelsize*Nx);
    Ini_PixelShifty(nL)=round(Ini_NAy(nL)*Pixelsize*Ny);
    
    v_s=Ini_NAx(nL);
    u_s=Ini_NAy(nL);

    % Generate shifted Green's function, spatial frequency terms
    G = real(1./(k_Medium*sqrt(1-lambda^2.*((uu-u_s).^2+(vv-v_s).^2))));
    Gf = real(1./(k_Medium*sqrt(1-lambda^2.*((uu+u_s).^2+(vv+v_s).^2))));
    uv_vector1=real(sqrt(1-lambda^2.*((uu-u_s).^2+(vv-v_s).^2)));
    uv_vector2=real(sqrt(1-lambda^2.*((uu+u_s).^2+(vv+v_s).^2)));
    uv_s=sqrt(1-lambda^2.*(u_s^2+v_s^2));

    % Shift pupil function to match illumination angle, complex conjugate
    Pupil = circshift(Aperture_fun,[Ini_PixelShiftx(nL),Ini_PixelShifty(nL)]);
    Pupilf = circshift(Aperture_fun,-[Ini_PixelShiftx(nL),Ini_PixelShifty(nL)]);
    
    for j=1:length(Depth_Set)
        tmp_PTF =...
            (Pupil  .* sin(k_Medium.*Depth_Set(j).*(uv_vector1 - uv_s)).*G+...
             Pupilf .* sin(k_Medium.*Depth_Set(j).*(uv_vector2 - uv_s)).*Gf)+...
         1i*(Pupil  .* cos(k_Medium.*Depth_Set(j).*(uv_vector1 - uv_s)).*G-...
             Pupilf .* cos(k_Medium.*Depth_Set(j).*(uv_vector2 - uv_s)).*Gf);
    
        tmp_ATF =...
            -(Pupil  .* cos(k_Medium.*Depth_Set(j).*(uv_vector1 - uv_s)).*G+...
              Pupilf .* cos(k_Medium.*Depth_Set(j).*(uv_vector2 - uv_s)).*Gf)+...
          1i*(Pupil  .* sin(k_Medium.*Depth_Set(j).*(uv_vector1 - uv_s)).*G-...
              Pupilf .* sin(k_Medium.*Depth_Set(j).*(uv_vector2 - uv_s)).*Gf);         

        PTF_4D(:, :, j, nL) = 0.5*dz*k_Wavenum^2.*fftshift(tmp_PTF);
        ATF_4D(:, :, j, nL) = 0.5*dz*k_Wavenum^2.*fftshift(tmp_ATF);
    end

end
clear tmp_PTF tmp_ATF

%% Invert scattering to recover the 3D Object
disp('Performing reconstruction...');
tic
for q = 1:2
    sum_PTF = 0;
    sum_ATF = 0;

    conj_PTF_Iten=0;
    conj_ATF_Iten=0;

    conj_term1=0;
    conj_term2=0;
    
    % Sort images into "Hot" and "Cold" images
    if(q == 1)
        I_tmp = I_Raw(:,:,1:2:(2*Length_MN));
    else
        I_tmp = I_Raw(:,:,2:2:(2*Length_MN));
    end

    for nL=1:Length_MN
        Itmp = I_tmp(:,:,nL);        
        Ihat_tmp = fft2((Itmp));
        
        sum_PTF=sum_PTF+abs(PTF_4D(:,:,:,nL)).^2;
        sum_ATF=sum_ATF+abs(ATF_4D(:,:,:,nL)).^2;

        conj_PTF_Iten=conj_PTF_Iten+conj(PTF_4D(:,:,:,nL)).*repmat(Ihat_tmp,1,1,length(Depth_Set));
        conj_ATF_Iten=conj_ATF_Iten+conj(ATF_4D(:,:,:,nL)).*repmat(Ihat_tmp,1,1,length(Depth_Set));

        conj_term1=conj_term1+conj(PTF_4D(:,:,:,nL)).*ATF_4D(:,:,:,nL);
        conj_term2=conj_term2+conj(ATF_4D(:,:,:,nL)).*PTF_4D(:,:,:,nL);
    end


    %% Apply inversion for Hot and Cold image sets
    if(q == 1)
        Normalized_term=(sum_PTF+Tau(1)).*(sum_ATF+Tau(2))-conj_term1.*conj_term2;
        v_re_hot = ((sum_ATF+Tau(2)).* conj_PTF_Iten - conj_term1.* conj_ATF_Iten) ./ Normalized_term;
        v_im_hot = ((sum_PTF+Tau(1)).* conj_ATF_Iten - conj_term2.* conj_PTF_Iten) ./ Normalized_term;
    else
        Normalized_term=(sum_PTF+Tau(1)).*(sum_ATF+Tau(2))-conj_term1.*conj_term2;
        v_re_cold = ((sum_ATF+Tau(2)).* conj_PTF_Iten - conj_term1.* conj_ATF_Iten) ./ Normalized_term;
        v_im_cold = ((sum_PTF+Tau(1)).* conj_ATF_Iten - conj_term2.* conj_PTF_Iten) ./ Normalized_term;
    end
end
toc

% Convert recovered permittivity values to the real space
v_im_hot = real(ifft2((v_im_hot)));
v_re_hot = real(ifft2((v_re_hot)));
v_im_cold = real(ifft2((v_im_cold)));
v_re_cold = real(ifft2((v_re_cold)));

% Convert complex permittivity to complex refractive index
n_re_hot = sqrt(((n_Medium.^2 + v_re_hot) + sqrt((n_Medium^2 + v_re_hot).^2 + v_im_hot.^2)) / 2);
n_im_hot = v_im_hot ./ n_re_hot ./ 2;
n_re_cold = sqrt(((n_Medium.^2 + v_re_cold) + sqrt((n_Medium^2 + v_re_cold).^2 + v_im_cold.^2)) / 2);
n_im_cold = v_im_cold ./ n_re_cold ./ 2;

RI_hot=n_re_hot+1i*n_im_hot;
RI_cold = n_re_cold + 1i*n_im_cold;

clear conj_term1 conj_term2 
clear n_re_hot n_re_cold n_im_hot n_im_cold v_re_hot v_im_hot v_re_cold v_im_cold sum_PTF sum_ATF conj_term1 conj_term2



