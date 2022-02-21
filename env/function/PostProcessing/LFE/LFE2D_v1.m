function [MU, PHI, LAMBDA, RR, CR] = LFE2D_v1(disp, dx, fms, varargin)
% Function LFE2D is to estimate G value based on local frequency estimate. No curl is applied
% Note: the NF is the number of filters applied; B is the filter bandwith;
% These two internal paramters are crucial for good G estimate.
%
% Input Variables:
% disp - displacement volume data with respect to time[row, col, timeFrames]
%   in microns. In the gel case, the displacement is along FH direction,
%   which is defined as z-direction, along the magnet bore axis.
% dx - pixel size in mm
% fms - vibration freqency in Hz
%
% Output Variables:
% MU - G* value
% PHI - phase angle
% LAMBDA - local wave length
% RR - local freuqency
% CR - confidence map
%
% Sep-01-2017===YF===Original Code
% Sep-02-2017===YF===Modify to a function
% Sep-06-2017===YF===Add optional code: certainty threshold and estimate output
% Sep-07-2017===YF===Add NF optional input

% addpath('G:\Postdoc research\MRE\LFE\updated_3Dcode');
% addpath('/Volumes/Suda/Postdoc research/MRE/LFE/updated_3Dcode')
% addpath('/Volumes/Suda/Suda Research/MRE/imProcess/LFE')
% addpath(genpath('/Volumes/Suda/Suda Research/Matlab/suda new'))

if nargin == 3
    % number of filters applied, smaller NF will result in a larger CR, but
    % with larger G estimates, larger NF will result in a smaller CR, but with
    % smaller G. NF = 11 is a tested number for reasonable estimates.
    NF = 15;
elseif nargin == 4
    NF = varargin{1};
end


plotfield = 0; % plot FFT principal component

% % sample data
% load('G:\Suda Research\MRE\experiment\MRE20170416\Phantom\DICOM\gelMRE20170416.mat');
% load('/Volumes/Suda/Suda Research/MRE/experiment/MRE20170416/Phantom/DICOM/gelMRE20170416.mat')

% displacement
Uz = disp*1e-3; % mm

%% FFT to get the 1st principal component
Fuz = fft(Uz, [], 3);
Uz1 = squeeze(Fuz(:,:,2));
% get real and imaginary components
% userows = 1:288; usecols = 1:288;
Wr=real(Uz1);
Wi=imag(Uz1);
% average filter smoothing, with a filter window size of 3x3
Crz = smooth2a(Wr,1,1) ;
Ciz = smooth2a(Wi,1,1);

% get rid of NaNs for LFE if they exist
Crz(find(isnan(Crz)))=zeros(size(find(isnan(Crz))));
Ciz(find(isnan(Ciz)))=zeros(size(find(isnan(Ciz))));

%% use complex field for inversion
LFARR = Crz;
LFARI = Ciz;
Ii = sqrt(-1);
LFARRtot = LFARR+Ii*LFARI;

if plotfield
    figure
    clim=max([max(abs(LFARR(:))),max(abs(LFARI(:)))]);
    subplot(1,2,1);
    temp = LFARR;
    imagesc(temp(:,:,1));
    axis image;axis off;colorbar;title('Real part')%caxis(ceil(10*clim)/10.*[-1 1]);
    subplot(1,2,2);
    temp = LFARI;
    imagesc(temp(:,:,1));
    axis image;axis off;colorbar;title('Imaginay part')
end

Xr = LFARRtot;
[RR,CR]=run_LFE_2d_v1(Xr, NF, 0);

%% POSTPROCESSING AND PLOTTING
if plotfield
    figure
    t1 = sprintf('|k| S%d',1);
    subplot(1,2,1)
    imagesc(RR,[min(RR(:)),max(RR(:))]); axis('off'); axis image; title(t1);
    t2 = sprintf('C S%d',1);
    subplot(1,2,2)
    imagesc(CR,[0 1]);axis('off');title(t2); axis image;
end;

Npix = size(disp,1);
voxsize = dx*1e-3; % m in the plane
L = Npix*voxsize;
FREQ = fms;
RHO = 1000; % kg/m^3
% convert local frequency to local wavelength
LAMBDA = L./(RR);       % using real component
%LAMBDA = L./(RI);       % using imag component
C = LAMBDA*FREQ; % m/s
MU = RHO*C.^2; % Pa

% LOSSINESS - ?
% estimate phase shift (loss)
LapLFARRtot = del2(LFARRtot);
PHI = angle(-LFARRtot./LapLFARRtot);

% optional code for output
% figure;
% subplot(1,3,1); imagesc(CR, [0 1]); axis image; colorbar; title('certainty')
% certaintyThreshold = 0.05;
% cMask = (CR.*mask)>certaintyThreshold;
% subplot(1,3,2); imagesc(cMask); axis image; title('certainty>0.05')
% subplot(1,3,3); imagesc(MU.*double(cMask)); axis image; colorbar; title('\mu')
%
% fprintf('mean MU = %f (Pa) +- %f (Pa)', mean(MU(cMask)), std(MU(cMask)));
