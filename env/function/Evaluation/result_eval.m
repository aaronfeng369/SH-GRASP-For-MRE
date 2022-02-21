function [psnr_pha,psnr_stiff,ssim_pha,ssim_stiff] = result_eval(img_pha,img_stiff,ref_pha,ref_stiff)
% ------------------------------------------------------
% Usage:
% used to evaluate the quality of reconstruction result
%
% Input: 
% (1) phase image after reconstruction and post-processing
% (2) stiffness image after reconstruction
% (3) reference phase image
% (4) reference stiffness image
% 
% Output:
% (1) psnr of phase and stiffness images
% (2) ssim of phase and stiffness images
% ------------------------------------------------------
% Runke Wang
% 2020/10/06 --- Original code
% 2022/01/20 --- Modified for stiffness evaluation
% ------------------------------------------------------
[x,y,z] = size(img_stiff);

[~,psnr_pha] = rmse_calc(img_pha,ref_pha,0);
[~,psnr_stiff] = rmse_calc(img_stiff,ref_stiff,0);
[~,ssim_pha_2d,~] = SSIM(img_pha,ref_pha);
[~,ssim_stiff_2d,~] = SSIM(img_stiff,ref_stiff);

% if z>1
%     ssim_pha = ssim_pha_3d;
%     ssim_mag = ssim_mag_3d;
% else
%     ssim_pha = ssim_pha_2d;
%     ssim_mag = ssim_mag_2d;
ssim_pha = ssim_pha_2d;
ssim_stiff = ssim_stiff_2d;
end