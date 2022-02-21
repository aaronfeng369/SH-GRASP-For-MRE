function [PNR] = PNR_cal_v2(Im_pos,Im_neg,mask_inp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-to-noise ratio (PNR) calculation for magnitude image
%
% Input Variables:
% Im_pos - Image encoded by MEG with positive polarity [row,column] (Complex value)
% Im_neg - Image encoded by MEG with negative polarity [row,column] (Complex value)
%
% Mask - mask for ignoring background
% 
% Output Variables:
%
% PNR - signal to noise ratio after correction 
%
% Reference:
% Bernstein MA, Ikezaki Y. Comparison of phase-difference and complex-difference processing in phase-contrast MR angiography. 
% J Magn Reson Imaging 1991;1(6):725-729
% 
% Record of Revisions:
% Apr-6-2021===RKW===Original Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% mask prepare
if nargin == 3
    se = strel('disk',5);
    [x,y,n_mask]=size(mask_inp);
    for n=1:n_mask
        mask_temp = ones(size(mask_inp(:,:,1)));
        mask_temp(isnan(mask_inp(:,:,n))) = 0;
        mask_temp = imerode(mask_temp, se);
        mask_temp(mask_temp==0) = NaN;
        mask(:,:,n) = mask_temp;
    end
else 
    n_mask=1;
    mask = ones(size(Im_pos));
end

%%
M_pos = abs(Im_pos);
M_neg = abs(Im_neg);

theta_pos = angle(Im_pos);
theta_neg = angle(Im_neg);

for n=1:n_mask
    A = (M_pos+M_neg)./2.*mask(:,:,n);
    del_the = unwrap_phase((theta_pos-theta_neg).*mask(:,:,n))./2;
    
    M(:,:,1)=M_pos;M(:,:,2)=M_neg;
    sig_M = nanstd(M,0,3).*mask(:,:,n);

    S = abs(A.*del_the.*2);
    N = ( 2* (sig_M.^2) .* (1+(del_the).^2) ).^(0.5);
    PNR_temp(:,:,n) = S./N;
    
end

PNR = nanmean(PNR_temp,'all');










