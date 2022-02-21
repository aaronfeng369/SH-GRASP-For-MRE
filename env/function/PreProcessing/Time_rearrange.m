function [data_re] = Time_rearrange(data_pos,data_nag,mode)
% ----------------------------------
% This function is used for rearranging the k spokes for external
% acquisition and generate the right k-space trajectory of wach phase
%
% Input:
% (1) Rawdata encoded by positive MEG
% (2) Rawdata encoded by nagative MEG
% (3) mode: 0 - concatenate directly; 
%           1 - rearrange to a form of continuous wave
%
% Output:
% Rearranged data
% ----------------------------------
% Runke Wang
% Version 2021/03/12
%-----------------------------------
n_dim = ndims(data_pos);
if mode == 0
    data_re=cat(n_dim,data_pos,data_nag);
elseif mode == 1
    data_nag_re = circshift(data_nag,2,n_dim);
    data_re = cat(n_dim,data_pos,data_nag_re);
    
end


end