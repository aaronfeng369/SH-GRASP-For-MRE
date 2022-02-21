function phase_filtered = DENSE_filter(result_phase)
% ------------------------------------------------------
% Usage:
% used for batch normalization after each external iteration
% Input: 
% (1) the result of reconstructed phase
% Output:
% (1) phase after filtering
% ------------------------------------------------------
% Record of Revisions
% Oct-19-2020===Runke Wang===Original Code
% Apr-10-2021===Yuan Feng===Minor modifications
% ------------------------------------------------------

[NPe,NFe,NFrame] = size(result_phase);
%gausian filter
phase_temp = zeros(NPe,NFe,1,NFrame);
phase_temp(:,:,1,:) = result_phase;
phase_temp = Filter3D(phase_temp);
phase_filtered = squeeze(phase_temp(:,:,1,:));