function [Gp, Gpp, nerr, lr2f ,li2f, r2, i2, r24, i24] = ...
    DI(dx,dy, NXY, mask, fms, NRO, NPH, filter, phi_array)
% 2-D linear equation fit, including complex constant (DC offset)
% spatially filter displacement and laplacian
% eliminate edges of ROIs (change imerode parameters to include more/less)
% 
% Input Variables:
% dx - pixel size(m)
% NXY - fitting area size, the area is +-NXY of the center pixel for fitting.
% mask - image mask
% fms - wave frequency(Hz)
% NRO - image readout number
% slice - image slice for calculation
% filter - .type1: ideal filter - 'i'; butterworth - 'b'
%          .band: D0 - passband
%          .type2: 'l'-lowpass, 'h'-high pass, 'b'-bandpass 
% 
% Output Variables:
% Gp - G'(Pa)
% Gpp - G''(Pa)
% nerr - normalized root meas square error 
% r24,i24 - real and imaginary parts of the FFT phase images after
% filtering.(M*N*NP)
% 
% Example: dx = .25e-3; fms = 400; NRO = 192;
% 
% Record of Revision
% Sep.22,2010=======EHK,PVB=======Original Code
% Jan.25,2010=======YAF=======Add r24, i24 output
% Feb.9 ,2010=======YAF=======Add NXY input
% Feb-21-2017===YF===modified for anisotropic voxels
% Mar-19-2017===YF===Add output display
% Jul-20-2017===YF===Correct the nan data that filter2 cannot deal with

% phase difference image
PDiff = squeeze(phi_array(:,:,:));
% angluar frequency 
w = 2*pi*fms;  
% get rid of ~2 edge voxels (erode mask)
% msk = imerode(mask,ones(5,5));  

% Initializing output
Gp = zeros(size(PDiff(:,:,1))); 
Gpp = Gp;
nerr = Gp;

% spatial filter, for ferret 'b', D0 = [0.01 0.99], for 1xgel 'i', D0 = [0.03 0.33]
if filter.type1 == 'i'
    [H, h, g] = ideal_filter_v2(NRO, NPH, filter.band, filter.type2); 
elseif filter.type1 == 'b'
    [H, h, g] = butterworth_filter_v2(NRO, NPH, filter.band, filter.order, filter.type2);
else
    error('please input a valid filter type, default: ideal')
    [H, h, g] = ideal_filter_v2(NRO, NPH, filter.band, filter.type2);  
end

% get Fourier coefficients
fphi = fft(PDiff,[],3);

% for 4 components
for ii = 1:size(fphi,3) 
    rphi(:,:,ii) = real(fphi(:,:,ii));
    iphi(:,:,ii) = imag(fphi(:,:,ii));
    
    % filter2 cannot deal with nans, must pad with zeros
    rphi(isnan(rphi)) = 0;
    iphi(isnan(iphi)) = 0;
     
    % filter displacement field
%     r24(:,:,ii) = filter2(h, rphi(:,:,ii).*mask);
%     i24(:,:,ii) = filter2(h, iphi(:,:,ii).*mask);
    r24(:,:,ii) = filter2(h, rphi(:,:,ii));
    i24(:,:,ii) = filter2(h, iphi(:,:,ii));
end
%figure
%mesh(rphi(:,:,2))
% get fundamental comoponents
r2 = r24(:,:,2); i2 = i24(:,:,2);
% the result is the same if we use the conjugate of 4th component
% r2 = r24(:,:,4); i2 = -i24(:,:,4);

% =========================inverse transform===============================
% % if we do not use the FFT component for the following calculation, we can
% % still inverse transform to get u in time domain for G calculation.
% uf2 = r2 + i*i2; 
% uf4 = r24(:,:,4)-i*i24(:,:,4); % we can not simply use conj for geting uf4 but do it in parts.
% uf_pad = zeros(size(PDiff)); 
% uf_pad(:,:,2) = uf2; 
% uf_pad(:,:,4) = uf4;
% u_filter = ifft(uf_pad,[],3);
% % the IFFT result has a phase shift so 1st and 3rd components can be used,
% % 2nd and 4th components are zero.Because we are using the u(x,t) for
% % spacial operation "del", so each time component are the same for
% % calculation, t is not an active parameter.
% r2 = real(u_filter(:,:,3)); i2 = imag(u_filter(:,:,3));
%==========================================================================


%=========================laplacian field==================================
% lr2f = 4*del2(r2,dx,dx);
% li2f = 4*del2(i2,dx,dx);
% we may use 3D laplacian field
% we may apply mask to laplacian field
lr2f = 4*del2(r2,dx,dy).*mask;
li2f = 4*del2(i2,dx,dy).*mask;
%==========================================================================

% figure(1); imagesc(msk>0); axis image; title('mask used')

% estimate modulus
% Prepare variables for fit to 2-D linear equation
density = 1000; % kg/m^3
% [II,JJ] = find(msk);
[II,JJ] = find((~isnan(mask)));


for ii=1:length(II),
    I0 = (II(ii)-NXY):(II(ii)+NXY);  %选取计算区域范围，NXY决定大小
    J0 = (JJ(ii)-NXY):(JJ(ii)+NXY);
    R2 = r2(I0,J0); %选取u在该范围内的数据
    I2 = i2(I0,J0);
    LR2 = lr2f(I0,J0);%del2(u)在该范围内的数据
    LI2 = li2f(I0,J0);

%=====================================
% % use square matrix equations
%     I0 = [II(ii), II(ii)+1];
%     J0 = [JJ(ii), JJ(ii)+1];
%     R2 = r2(I0,J0);
%     I2 = i2(I0,J0);
%     LR2 = lr2f(I0,J0);
%     LI2 = li2f(I0,J0);
%=====================================
    x=-density*w^2*R2(:);  %排列成列向量
    y=-density*w^2*I2(:);
    a=LR2(:);
    b=LI2(:);

    %==============================
    % why multiply by std(b)
%     u=ones(size(b))*std(b);
    u=ones(size(b));
     %==============================
     
    z=zeros(size(b));


    % Algebra
    % x + i*y = (g1 + i*g2)*(a + i*b) + (c1 + i*c2);
    % Real and imaginary parts
    % x = g1*a - g2*b + c1;
    % y = g1*b + g2*a + c2;

    % Perform the fit via pseudoinverse
    xy=[x;y];           % 2Nx1
    ab=[a,-b,u,z;
        b,a,z,u];      % 2Nx4    (xy = ab*G)   ( G=[g1;g2;c1;c2] )
    G = ab\xy;

    % Check - fitted model is
    XY = ab*G;
    xch2=XY(1:length(b));
    ych2=XY(length(b)+1:end);


    % complex modulus and constants
    Gp(II(ii),JJ(ii)) = G(1);
    Gpp(II(ii),JJ(ii)) = G(2);

    %normalized rms err
    nerr(II(ii),JJ(ii)) = sqrt(sum((ych2-y).^2 + (xch2-x).^2)) / sqrt(sum((y).^2 + (x).^2));

end
Gp(isnan(Gp)) = 0;
Gpp(isnan(Gpp)) = 0;
% indMask = (nerr<0.5) & (Gp >0) & (Gpp>0) & msk;
indMask = (nerr<0.5) & (Gp >0) & (Gpp>0) & ~isnan(mask);
% 
ok = find(indMask);
% 
% figure
% subplot(2,2,1); imagesc(Gp); colormap('jet'); axis image; colorbar; 
% title(['G''(Pa), fms =', num2str(fms), 'Hz'], 'fontsize', 22, 'fontname', 'Garamond', 'FontWeight', 'bold')
% 
% subplot(2,2,2); imagesc(Gpp);colormap('jet'); axis image; colorbar; 
% title(['G''''(Pa), fms =', num2str(fms), 'Hz'], 'fontsize', 22, 'fontname', 'Garamond', 'FontWeight', 'bold')
% 
% subplot(2,2,3); imagesc(nerr,[0 1]); colormap('jet'); axis image; colorbar;
% title('normalized RMS error', 'fontsize', 22, 'fontname', 'Garamond', 'FontWeight', 'bold')
% 
% G_star = sqrt((nanmean((Gp(ok))))^2+(nanmean((Gpp(ok))))^2);
% fprintf('G'' = %f +- %f Pa\n', nanmean((Gp(ok))), nanstd(Gp(ok)));
% fprintf('G'''' = %f +- %f Pa\n', nanmean((Gpp(ok))), nanstd(Gpp(ok)));
% fprintf('G* = %f Pa\n', G_star);
% 
% figure; 
% imagesc(indMask); axis image; colormap gray;
% title('(nerr<0.5) & (Gp >0) & (Gpp>0) & msk')

end