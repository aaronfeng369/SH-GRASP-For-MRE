function [rawdata_re,k_traj_re] = k_rearrange_v3(rawdata,param_nufft)
% ----------------------------------
% Runke Wang
% Version 2021/01/16
% Version 2021/11/13 : Add signal intensity calibration
% ----------------------------------
% This function is used for rearranging the k spokes for external
% acquisition and generate the right k-space trajectory of wach phase
%
% Input:
% original k space 
% [NPe,NSPe,NAvg,Necho,NFe,Ncoil]=size(rawdata)
% Npe - spokes num
% NAvg - phase num
% NFe - sampling point num of readout
% Ncoil - coil num
%
% Output:
% k-spokes of each phase after rearrangement.
% operaterE for nufft
%-----------------------------------

% Signal Intensity Calibration
% rawdata = squeeze(rawdata);
% [NPe,NFe,Ncoil]=size(rawdata);
% for nch = 1:Ncoil
%     for npe = 1:NPe
%         rawdata(npe,:,nch) =  rawdata(npe,:,nch)/max(abs(rawdata(npe,:,nch)));
%     end
% end

% Generate k-traj for whole
rawdata=permute(squeeze(rawdata),[2,1,3]);
[NFe,NPe,Ncoil]=size(rawdata);
Nphase = param_nufft.Nphase;%nt=NPe/Nphase;
k_traj_whole = trajRadialGR(param_nufft.RO_point,param_nufft.Nspokes,param_nufft.os_rate);


% re-arrange for k-sapce and generate k-trajectory
nt = NPe/Nphase;
rawdata_re = zeros(NFe,nt,Ncoil,Nphase);
k_traj_re = zeros(NFe,nt,Nphase);
for p=1:Nphase
    for c=1:Ncoil
        rawdata_re(:,:,c,p) = rawdata(:,p:Nphase:end,c);
    end
    k_traj_re(:,:,p) = k_traj_whole(:,p:Nphase:end);
end


end