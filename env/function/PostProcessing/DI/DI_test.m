%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for figure 3 (MRD part):
% Reconstruction for 4 phase images of phantom acquired by MRD
%
% Usage:
% (1) Change the EnvPath to the path of the file 'env'
% (2) Change DataPath_read to the location of MRD data
%
% Outputs:
% (1) 4 reconstucted phase images acquired by MRD for part of figure 3.
% (2) |G*|map of the phantom acquired by MRD
% (3) mean mouluses of different components in the phantom
%
% Record of Revisions:
% March-09-2021===RKW===Original Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
%% -------------Add path---------- 
EnvPath = 'F:\new_algorithm\FC_sparse_SENSE\图片绘制\new\code\代码整理\env\'; % Please change this path while using
addpath ([EnvPath,'function\']);
addpath ([EnvPath,'FengYuanLab_toolbox\']);
addpath ([EnvPath,'FengYuanLab_toolbox\common\']);
addpath ([EnvPath,'FengYuanLab_toolbox\for_LFE\']);
addpath ([EnvPath,'grasp_v2\']);
addpath ([EnvPath,'grasp_v2\nufft_toolbox\']);
load([EnvPath,'cmapAWAVE.mat']);

K=2*pi*(42.58*10^6)*(40*10^(-3))*((250+800)*10^(-6)); % Set motion sensitivity
x_target = 128;y_target = 128; % Set the data size of output images

%% -------------Read data---------- 
DataPath_read = 'F:\new_algorithm\FC_sparse_SENSE\图片绘制\new\code\figure5_sequence_comp_phantom\data\20210117\phantom_20hz\'; % Please change this path while using
DataPath_save = [pwd,'\result_temp\'];

%load aquired radial data
filename_1 = [DataPath_read,'rawdata\UID_1610887223_MRD_v2b_180_pos.raw']; 
filename_2 = [DataPath_read,'rawdata\UID_1610887259_MRD_v2b_180_nag.raw'];
rawdata_1 = Read_UIH_Raw_v5_NAVG(filename_1); 
rawdata_2 = Read_UIH_Raw_v5_NAVG(filename_2);
[NPe,NSPe,NAvg,Necho,NFe,Ncoil]=size(rawdata_1)

% Load prior information (Sensitivity map and mask)
prior_infor = load([DataPath_read,'prior_infor\prior_infor_phantom_20hz.mat']);
Sensitivity_map = squeeze(prior_infor.Sensitivity_map(:,:,:,1));
Sensitivity_map = Sensitivity_map./max(abs(Sensitivity_map(:)));
maskNaN3D = prior_infor.maskNaN3D;
maskNaN3D_rot = rot90(crop(maskNaN3D,x_target,y_target),-1);

%% -------------Phase images reconstruction---------- 
% Initialization for nufft  
param_nufft.Nspokes = NPe;                        
param_nufft.RO_point = NFe/2;
param_nufft.os_rate = 2;
param_nufft.Nphase = 4;
NFrame = param_nufft.Nphase;

% Prepare k-space re-arrange and nufft operator
[rawdata_re_1,k_traj_re_1] = k_rearrange_v3(rawdata_1,param_nufft);
[rawdata_re_2,k_traj_re_2] = k_rearrange_v3(rawdata_2,param_nufft);

% Data concatenation for simultaneous reconstruction
kdata(:,:,:,1:4) = rawdata_re_1;
k_traj(:,:,1:4) = k_traj_re_1;
kdata(:,:,:,5:8) = rawdata_re_2;
k_traj(:,:,5:8) = k_traj_re_2;

% Density compensation
w_ref = abs(k_traj);
for phase=1:NFrame*2
    for ch=1:Ncoil
        kdata(:,:,ch,phase)=kdata(:,:,ch,phase).*sqrt(w_ref(:,:,phase));
    end
end

% Images reconstruction
operaterE=MCNUFFT(k_traj,abs(k_traj),Sensitivity_map);
recon_nufft_fs = operaterE' * kdata;

% Normalzation
recon_nufft_fs = recon_nufft_fs./max(recon_nufft_fs(:));

% Positive and negative cancelling, and phase unwrapping
phase_ref = zeros(NFe,NFe,NFrame);
for phase = 1:NFrame 
    phase_ref(:,:,phase) = unwrap_phase(angle(maskNaN3D.*((recon_nufft_fs(:,:,phase)./recon_nufft_fs(:,:,phase+4)))));
end

% Post-rocessing and display
[ref,ref_cat] = post_processing(phase_ref,recon_nufft_fs,K,[DataPath_save,'phase'],AWAVE);

%% Elastic modulus estimation (LFE)
% Calculate elastic modulus
[MU, PHI, LAMBDA, RR, CR] = LFE2D_v1(ref./(1*10^(-6)), 2.34, 20);
% Display |G*| map
figure(10);
imshow(maskNaN3D_rot.*MU,[],'border','tight','initialmagnification','fit');set(gcf,'Position',[0,0,x_target,y_target]);colormap jet;
caxis([0,7000]);
saveas(10, [DataPath_save,'G.bmp']);

%% Elastic modulus estimation (DI)
dx = 2.34; dy = 2.34; fms = 40;
% DI inversion
filter.type1 = 'i';
filter.type2 = 'b';
filter.band = [0.01 0.99]; % [0.03 0.3]; %
filter.order = 1;
NXY = 0;
NRO = size(ref,1);
NPH = size(ref,2); 
for ss = 1:1
    %dispt = squeeze(phase(:,:,ss,:,3,1)).*CF;
    %maskNaN = maskNaN_3d(:,:,ss);
    % DI inversion
    [Gp, Gpp, nerr, lr2f ,li2f, r2, i2, r24, i24] = ...
        DI(dx, dx, NXY, maskNaN3D_rot, fms, NRO, NPH, filter, ref);
    close all
    G_star(:,:,ss) = sqrt(Gp.^2.+Gpp.^2);
end

figure(11);
imshow(G_star,[],'border','tight','initialmagnification','fit');set(gcf,'Position',[0,0,x_target,y_target]);colormap jet;

