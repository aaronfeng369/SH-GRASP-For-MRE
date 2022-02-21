function mean_phase = meanphase_temporal(temp_image)
% ----------------------------------
% This function is used for calculate the mean phase for the phase image of MRE 
% Input:
% Image with DC component through time axis (3d data: x,y,t)
% Output:
% temporal mean phase of the ROI
% 
% Record of Revisions
% Oct-12-2020===Runke Wang===Original Code
% Apr-10-2021===Yuan Feng===Minor modifications
%-----------------------------------
temp_image = squeeze(temp_image);
[NFe,NPe,Nt] = size(temp_image);
phase_offset = zeros(NFe,NPe);
mean_phase = nanmean(temp_image,3);
% phase_offset = sum(temp_image,3);
% mean_phase = phase_offset./Nt;
