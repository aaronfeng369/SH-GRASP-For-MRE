%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for figure 5 : Reconstruction for simulation data
% 
% Usage:
% (1) Change the EnvPath to the path of the file 'env'
% (2) Set the reduction factor you want to use
%
% Outputs:
% (1) 4 phase and 8 magnitude images reconstructed by NUFFT
% (2) 4 phase and 8 magnitude images reconstructed by GRASP
% (2) 4 phase and 8 magnitude images reconstructed by CS-wavelet
% (3) 4 phase and 8 magnitude images reconstructed by SH-GRASP
% (4) RMSE and SSIM for evaluation
%
% Record of Revisions:
% Mar-09-2021===RKW===Original Code
% Jan-20-2022===RKW===Add CS-wavelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;

%% -------------Add path---------- 
EnvPath = 'F:\new_algorithm\FC_sparse_SENSE\图片绘制\new\code\代码整理\Github\env\'; % Please change this path while using
addpath (genpath(EnvPath));
load([EnvPath,'cmapAWAVE.mat']);

R = 20; % Please change the reduction factor 

CF = 1e-3*2*20/(42.58*15*2*pi);% m/rad,phase-to-displacement ratio 
K=1/(CF*4);
%% -------------Read data---------- 
DataPath_read = [EnvPath,'\data\']; 
DataPath_save = [pwd,'\result_temp\'];

% Load aquired radial data
load([DataPath_read,'brain_20hz.mat']); 
[Nfe,Npe,Ncoil,Nphase]=size(kdata_rad_pos);

% Load mask and sensitivity map
load([DataPath_read,'prior_infor.mat']);
maskNaN3D_rot = rot90(maskNaN3D,-1);

[Nx,Ny,Nc]=size(Sensitivity_map)
%% data reconstruction for reference
Current_path = 'ref';

w1 = abs(k_traj_pos);
w2 = abs(k_traj_nag);
for phase=1:Nphase
    for ch=1:Ncoil
        kdata_rad_pos(:,:,ch,phase)=kdata_rad_pos(:,:,ch,phase).*sqrt(w1(:,:,phase));
        kdata_rad_nag(:,:,ch,phase)=kdata_rad_nag(:,:,ch,phase).*sqrt(w2(:,:,phase));
    end
end

% Fully-sampled data reconstruction
operaterE_pos = MCNUFFT(k_traj_pos,abs(k_traj_pos),Sensitivity_map);  
im_ref_pos = operaterE_pos' * kdata_rad_pos;
operaterE_nag = MCNUFFT(k_traj_nag,abs(k_traj_nag),Sensitivity_map);  
im_ref_nag = operaterE_nag' * kdata_rad_nag;

% Normalzation in image domain
im_ref_pos = im_ref_pos./max([im_ref_pos(:);im_ref_nag(:)]);
im_ref_nag = im_ref_nag./max([im_ref_pos(:);im_ref_nag(:)]);
kdata_pos = operaterE_pos * im_ref_pos;
kdata_nag = operaterE_nag * im_ref_nag;

% Phase difference reconstruction
im_ref_pha = zeros(Nx,Ny,Nphase);
for phase = 1:Nphase 
    im_ref_pha(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((im_ref_pos(:,:,phase)./im_ref_nag(:,:,phase)))));
end

% post-rocessing and display
[img_fs_pha,~] = post_processing(im_ref_pha,im_ref_pos,K,[DataPath_save,Current_path,'\ref_'],AWAVE);

[Stiffness, PHI, LAMBDA, RR, CR] = LFE2D_v1(img_fs_pha./(1*10^(-6)), 2.34, 20);
filter = fspecial('average', 3);
Stiffness_ref = imfilter(maskNaN3D_rot.*Stiffness,filter);
figure(101);
imshow(Stiffness_ref,'border','tight','initialmagnification','fit');colormap hot;set(gcf,'Position',[0,0,128,128]);
caxis([0,3500]);
saveas(101, [DataPath_save,Current_path,'\Stiffness_ref','.bmp']);
%% retrospectively down-sampled data reconstruction and normalization
Current_path = 'ds';

Npe_ds = ceil(Npe/R); 
% Data down-sampling
kdata_pos_ds = kdata_pos(:,1:Npe_ds,:,:);
k_traj_pos_ds = k_traj_pos(:,1:Npe_ds,:);
kdata_nag_ds = kdata_nag(:,1:Npe_ds,:,:);
k_traj_nag_ds = k_traj_nag(:,1:Npe_ds,:);

% Down-sampled data reconstruction
operaterE_pos_ds = MCNUFFT(k_traj_pos_ds,abs(k_traj_pos_ds),Sensitivity_map);  
img_pos_ds = operaterE_pos_ds' * kdata_pos_ds;
operaterE_nag_ds = MCNUFFT(k_traj_nag_ds,abs(k_traj_nag_ds),Sensitivity_map);  
img_nag_ds = operaterE_nag_ds' * kdata_nag_ds;

% Normalzation in image domain
img_pos_ds = img_pos_ds./max([img_pos_ds(:);img_nag_ds(:)]);
img_nag_ds = img_nag_ds./max([img_pos_ds(:);img_nag_ds(:)]);


im_ds_pha = zeros(Nx,Ny,Nphase);
for phase = 1:Nphase
    im_ds_pha(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((img_pos_ds(:,:,phase)./img_nag_ds(:,:,phase)))));
end


% post-rocessing and display
[img_ds_pha,~] = post_processing(im_ds_pha,img_pos_ds,K,[DataPath_save,Current_path,'\ds_R',num2str(R),'_'],AWAVE);
img_ds_pha = img_ds_pha.*K.*2; img_ref_pha = img_fs_pha.*K.*2; % Compared in phase domain after post-processing;

[Stiffness, PHI, LAMBDA, RR, CR] = LFE2D_v1(img_ds_pha./(1*10^(-6)), 2.34, 20);
filter = fspecial('average', 3);
Stiffness_ds = imfilter(maskNaN3D_rot.*Stiffness,filter);
figure(201);
imshow(Stiffness_ds,'border','tight','initialmagnification','fit');colormap hot;set(gcf,'Position',[0,0,128,128]);
caxis([0,3500]);
saveas(201, [DataPath_save,Current_path,'\Stiffness_ds_R',num2str(R),'.bmp']);

[psnr_stiffness_zp,psnr_phase_zp,ssim_stiffness_zp,ssim_phase_zp] = result_eval(Stiffness_ds,img_ds_pha,Stiffness_ref,img_ref_pha);
%% Initialization for parameter of CG
param.lambda=0;
param.sTVWeight = 0;
param.lambda_fft=0;
param.lambda_mreTV_L1=0;
param.lambda_mreTV_L2=0;
param.Denoising_on = 0;
param.DC_removing_on = 0;
param.WaveletWeight = 0;
param.nite = 5;
param.display=1;

% Sparsifying transform
param.W = TV_Temp();
param.sTV = TVOP_v2();
param.wavelet = Wavelet('Daubechies',4,4);
%% ------------- iGRASP_tTV reconstruction ----------
Current_path = 'grasp';
% Parameter setting for GRASP
param.lambda=0.02;
param_pos = param;param_nag = param;
param_pos.E = operaterE_pos_ds;param_nag.E = operaterE_nag_ds;
param_pos.y = kdata_pos_ds;param_nag.y = kdata_nag_ds;

% reconstruction
fprintf('\n GRASP reconstruction -- tTV\n')
tic
im_grasp_pos=img_pos_ds;
im_grasp_nag=img_nag_ds;
for n=1:5
	im_grasp_pos = CSL1NlCg_mre_v4(im_grasp_pos,param_pos);
    im_grasp_nag = CSL1NlCg_mre_v4(im_grasp_nag,param_nag);
end
toc

im_grasp_pha = zeros(Nx,Ny,Nphase);
for phase = 1:Nphase 
    im_grasp_pha(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((im_grasp_pos(:,:,phase)./im_grasp_nag(:,:,phase)))));
end

% Post-rocessing and evaluation
[img_grasp_pha,~] = post_processing(im_grasp_pha,im_grasp_pos,K,[DataPath_save,Current_path,'\GRASP_R',num2str(R),'_'],AWAVE);
img_grasp_pha = img_grasp_pha.*K.*2;

[Stiffness, PHI, LAMBDA, RR, CR] = LFE2D_v1(img_grasp_pha./(1*10^(-6)), 2.34, 20);
filter = fspecial('average', 3);
Stiffness_grasp = imfilter(maskNaN3D_rot.*Stiffness,filter);
figure(301);
imshow(Stiffness_grasp,'border','tight','initialmagnification','fit');colormap hot;set(gcf,'Position',[0,0,128,128]);
caxis([0,3500]);
saveas(301, [DataPath_save,Current_path,'\Stiffness_grasp_R',num2str(R),'.bmp']);

[psnr_stiffness_grasp,psnr_phase_grasp,ssim_stiffness_grasp,ssim_phase_grasp] = result_eval(Stiffness_grasp,img_grasp_pha,Stiffness_ref,img_ref_pha);
%% CS-wavelet
Current_path = 'cs';
% Parameter setting for CS
param.lambda=0.0;
param.WaveletWeight = 0.02;
param_pos = param;param_nag = param;
param_pos.E = operaterE_pos_ds;param_nag.E = operaterE_nag_ds;
param_pos.y = kdata_pos_ds;param_nag.y = kdata_nag_ds;

% reconstruction
fprintf('\n cs_wavelet reconstruction \n')
tic
im_cs_pos=img_pos_ds;
im_cs_nag=img_nag_ds;
for n=1:2
	im_cs_pos = CSL1NlCg_mre_v4(im_cs_pos,param_pos);
    im_cs_nag = CSL1NlCg_mre_v4(im_cs_nag,param_nag);
end
toc

im_grasp_pha = zeros(Nx,Ny,Nphase);
for phase = 1:Nphase 
    im_cs_pha(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((im_cs_pos(:,:,phase)./im_cs_nag(:,:,phase)))));
end

% Post-rocessing and evaluation
[img_cs_pha,~] = post_processing(im_cs_pha,im_cs_pos,K,[DataPath_save,Current_path,'\CS_R',num2str(R),'_'],AWAVE);
img_cs_pha = img_cs_pha.*K.*2;

[Stiffness, PHI, LAMBDA, RR, CR] = LFE2D_v1(img_cs_pha./(1*10^(-6)), 2.34, 20);
filter = fspecial('average', 3);
Stiffness_cs = imfilter(maskNaN3D_rot.*Stiffness,filter);
figure(501);
imshow(Stiffness_cs,'border','tight','initialmagnification','fit');colormap hot;set(gcf,'Position',[0,0,128,128]);
caxis([0,3500]);
saveas(301, [DataPath_save,Current_path,'\Stiffness_cs_r',num2str(R),'.bmp']);

[psnr_stiffness_cs,psnr_phase_cs,ssim_stiffness_cs,ssim_phase_cs] = result_eval(Stiffness_cs,img_cs_pha,Stiffness_ref,img_ref_pha);
%% ------------- proposed----------
Current_path = 'shgrasp';
% Rearrangement
im_shgrasp_0 = Time_rearrange(img_pos_ds,img_nag_ds,1);
kdata_ds = Time_rearrange(kdata_pos_ds,kdata_nag_ds,1);
k_traj_ds = Time_rearrange(k_traj_pos_ds,k_traj_nag_ds,1);
operaterE_ds = MCNUFFT(k_traj_ds,abs(k_traj_ds),Sensitivity_map);

% Use mre TV as a new regularization term
param.E = operaterE_ds;
param.y = kdata_ds;
param.lambda=0;
param.WaveletWeight=0;
param.lambda_mreTV_L1 = 0.02;
param.lambda_fft = 0.0001;

% reconstruction
fprintf('\n SH-GRASP reconstruction \n')
tic
im_shgrasp=im_shgrasp_0;
for n=1:5
	[im_shgrasp,cost(:,n)] = CSL1NlCg_mre_v4(im_shgrasp,param);
end
toc

im_shgrasp_pha = zeros(Nx,Ny,Nphase);
for phase = 1:Nphase 
    if phase <=2
        im_shgrasp_pha(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((im_shgrasp(:,:,phase)./im_shgrasp(:,:,phase+6)))));
    else
        im_shgrasp_pha(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((im_shgrasp(:,:,phase)./im_shgrasp(:,:,phase+2)))));
    end
end


% post-rocessing and display
[img_shgrasp_pha,~] = post_processing(im_shgrasp_pha,im_shgrasp,K,[DataPath_save,Current_path,'\SHGRASP_R',num2str(R),'_'],AWAVE);
img_shgrasp_pha = img_shgrasp_pha.*K.*2;

[Stiffness, PHI, LAMBDA, RR, CR] = LFE2D_v1(img_shgrasp_pha./(1*10^(-6)), 2.34, 20);
filter = fspecial('average', 3);
Stiffness_shgrasp = imfilter(maskNaN3D_rot.*Stiffness,filter);
figure(401);
imshow(Stiffness_shgrasp,'border','tight','initialmagnification','fit');colormap hot;set(gcf,'Position',[0,0,128,128]);
caxis([0,3500]);
saveas(401, [DataPath_save,Current_path,'\Stiffness_shgrasp_R',num2str(R),'.bmp']);

[psnr_stiffness_shgrasp,psnr_phase_shgrasp,ssim_stiffness_shgrasp,ssim_phase_shgrasp] = result_eval(Stiffness_shgrasp,img_shgrasp_pha,Stiffness_ref,img_ref_pha);