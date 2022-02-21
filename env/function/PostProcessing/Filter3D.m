function [PhiDCFilter] = Filter3D(PhiDC)
% The function file is to spatially filter data using smooth3: 
% 3D Gaussian spatial filtering
% 
% Input variables:
% PhiDC - DC corrected data after 3D unwrapping [PE RO slice# Phase#]
% 
% Output variables:
% The following arrays have the same structure : [PE RO slice# Phase#]
% PhiDCFilter - 3D Gaussian filtered data
% 
% Record of Revision
% Jan-18-2012===YAF===Original Code

%% 3D Gaussian filter 
% filter radius
[Npe,Nfe,Nslice,Nphase] = size(PhiDC);
nfpix = [1 2 3 4]; 
mm = 2; % original:1
% filter size
filter_sz = [2*nfpix(mm)+1 2*nfpix(mm)+1 2*nfpix(mm)+1]; 
% filter standard deviation
filter_std = nfpix(mm); 
% 3D Gaussian filter
PhiDCFilter = zeros(Npe,Nfe,2,Nphase);
for ii = 1:size(PhiDC,4)
    tic;
    phi_temp = zeros(Npe,Nfe,2);
    phi_temp(:,:,1) = PhiDC(:,:,:,ii);
    phi_temp(:,:,2) = PhiDC(:,:,:,ii);% keep the 3d dimension
    PhiDCFilter(:,:,:,ii) = smooth3(phi_temp,'gaussian',filter_sz,filter_std);
    disp(['Time to smooth: ',num2str(toc),'sec, filt_sz/std: ',num2str([filter_sz filter_std])]);
end
