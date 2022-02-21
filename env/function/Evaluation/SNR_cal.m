%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR calculation for magnitude image
%
% Input Variables:
% Im - magnitude image [row,column]
%
% Mask - mask for ignoring background
% 
% Output Variables:
% S - signal strength
%
% N - noise strength
%
% SNR - signal to noise ratio after correction 
%
% Record of Revisions:
% Nov-27-2020===SQ===Original Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,N,SNR] = SNR_cal(Im,Mask)
edge_ignore = 4;% ignore several pixels near the edge
NoisyRegionL = 20; % size of noise region
C = 0.7; %correction ratio for multi-coil
Mask01 = ones(size(Mask));
%Mask01(isnan(Mask01)) = 0; %original
Mask01(isnan(Mask)) = 0;
% erode mask
se = strel('disk',10);
maskerode = imerode(Mask01, se);
maskNaN_signal = ones(size(Mask));% signal region
maskNaN_noise = zeros(size(Mask));%nosie region
mask_dia = zeros(size(Mask));% diagram of all regions
% mask for SNR calculation
maskNaN_signal(maskerode==0) = NaN;
% four corner of the noise
maskNaN_noise(edge_ignore:edge_ignore+NoisyRegionL,edge_ignore:edge_ignore+NoisyRegionL) = 1;
maskNaN_noise(end-edge_ignore-NoisyRegionL:end-edge_ignore,edge_ignore:edge_ignore+NoisyRegionL) = 1;
maskNaN_noise(edge_ignore:edge_ignore+NoisyRegionL,end-edge_ignore-NoisyRegionL:end-edge_ignore) = 1;
maskNaN_noise(end-edge_ignore-NoisyRegionL:end-edge_ignore,end-edge_ignore-NoisyRegionL:end-edge_ignore) = 1;
maskNaN_noise(maskNaN_noise==0) = NaN;
mask_dia(~isnan(maskNaN_noise)) = 1;
%mask_dia(~isnan(maskNaN_signal)) = 1;

% figure
% imshow(mask_dia)
% figure
% imshow(Im,[])

% mean signal value
S_temp = nanmean(nanmean(Im.*maskNaN_signal)); % original

% std of noise value
N_temp1 = Im.*maskNaN_noise;
N_temp = nanstd(N_temp1(:)); % original
SNR_temp = S_temp/(N_temp/C); % original

SNR = SNR_temp;
S = S_temp;
N = N_temp;
end