function [PhiDCFilter] = FilterTemp(PhiDC, N2)
% The function file is to spatially filter data using smooth3: 
% 3D Gaussian spatial filtering
% 
% Input variables:
% PhiDC - DC corrected data after 3D unwrapping [PE RO slice# Phase#]
% N2 - number of time points in interpolated data (optional)
% 
% Output variables:
% The following arrays have the same structure : [PE RO slice# Phase#]
% PhiDCFilter - Temporally filtered data, only the 2nd and 4th component of
% FFTed input data were kept.
% 
% Record of Revision
% Jan-18-2012===YAF===Original Code

%% Temporal filter 

% fft to in the temporal domain to get fundamental spectrum
fphi = fft(PhiDC,[],4);

% number of time points in original data (input)
N1 = size(PhiDC,4);

% extract real and imaginary parts
rphi = real(fphi);      % real part (in-phase)
iphi = imag(fphi);      % imag part (out-of-phase)

% look at coefficients of fundamental component (in total )
% The spectrum of the FFT results are symmetrical, so the second(f1) and
% last (-f1) components are the "fundamental" componet, the first one is DC(f0=0).
% We only care about the first no-DC components now because the higher
% frequency ones may due to the noise or error that we may not control now.
r2  = rphi(:,:,:,2);
i2 =  iphi(:,:,:,2);
r4  = rphi(:,:,:,size(PhiDC,4));
i4 =  iphi(:,:,:,size(PhiDC,4));

% first initializing
fphi_filter = zeros(size(PhiDC));

% keep only the fundamental components
fphi_filter(:,:,:,2) = r2 + 1i*i2;
fphi_filter(:,:,:,N2) = r4 + 1i*i4;

% inverse FFT using fitlered components, since the points were multiplied
% by N2/N1, the IFT algorithm place the 1/M multiplier at the inverse part,
% so we need to compensate the multiplier by N2/N1
PhiDCFilter = (N2/N1)*ifft(fphi_filter, [], 4);
