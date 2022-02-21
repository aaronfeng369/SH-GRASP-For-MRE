function [xx,ff,Q,D,Ri] = locfreq_fine_2d(x1,ri,NK)
% Make harmonic phase image

% Input Variables:
% x1 = image matrix, r0 = filter radius, 
% NK =  direction (1, 2, or 3)
% 
% Outputs Variables:
% xx = filtered image, ff = filtered k-space, Q filter, 
% D = directional filter, Ri = radial filter
% 
% Record of Revisions
% line 39 changed by RJO on 8/9/13 and corrected 8/12/13
% Sep-01-2017===YF===Modified for 2D cases, change filter bandwith B to
%                       2*sqrt(2)

% Basic statistics
[NN,NC] = size(x1);
[ii,jj] = meshgrid(1:NC,1:NN);
cen = [NC/2+1;NN/2+1];
u = zeros(NN,NC,2);

% FFT and shifts data
f1 = fftn(x1);
f1(1) = 0; % DC peak
f1 = fftshift(f1);

% Define radial filter    (Eq 6)
% B = 2;  % Half octave bandwidth
B = 2*sqrt(2); % A center frequency ratio of 2, with a relative filter bandwidth of 2*sqrt(2)
Cb = 4/(B^2*log(2));
u(:,:,1) = ii-cen(1);
u(:,:,2) = jj-cen(2);
% u(:,:,:,3) = kk-cen(3);
rr = sqrt((ii-cen(1)).^2 + (jj-cen(2)).^2);
Ri = exp(-Cb*(log(rr/ri)).^2);

% directional filter - applied in x, y and z directions (Eqs.8-12)
D=zeros(size(Ri));
nk=zeros(2,1);
nk(NK)=1;
for m=1:NN
    for n=1:NC
        U = reshape(u(m,n,:),2,1);
        nU = norm(U);
        if nU==0,nU=1;end;
        D(m,n) = abs(U'*nk).*(norm(U'*nk)>0)/nU; %changed by RJO
    end;
end;

% total filter = radial* directional
Q = Ri.*D;

% filter
ff = Q.*f1;

% Inverse Fourier transform
xx=ifftn(ifftshift(ff));

xx=abs(xx); % why take absolute value here?
