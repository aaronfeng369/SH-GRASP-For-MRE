function [r,c,xxmat]=run_LFE_2d_v1(X, NF, debug)
% Run local filter repeatedly to get local frequency
% need image volume X
% Knutsson Westin Granlund IEEE 1994
% basic frequency (x 8 octaves)
% 
% Input Variables:
% X - Complex or real images of the principal component matrix from FFT
% NF - number of quadrature filters applied, a larger number of NF will get
%       a more confident estimate.; for 2D case, NF = 50 is good enough.
% 
% Output Variables:
% r - local frequency
% c - confidence map
% 
% Record of Revisions:
% modified by R. Okamoto 8/12/13 to call locfreq_fine_3d_rjoedit.m
% Sep-01-2017===YF===Modified for 2D cases

%mfilename  %display name of running function
r0=1;
[NR,NC]=size(X);
if nargin ==1, 
    debug=0;
end
if debug,
    xxmat = zeros(NR,NC,2,NF);
end

% wide range frequency estimate
for ii=1:NF,                                        % each iteration
    ri = (2^(ii/2))*r0;                             % center frequency
    q{ii}=zeros(size(X));
    for NK=1:2,                                     % three directions
        [xx,ff,Q,D,Ri] = locfreq_fine_2d(X,ri,NK);     % implement filter
        q{ii}=q{ii}+xx; 
        if debug, % chged by rjo to store intermediate results
            xxmat(:,:,:,NK,ii);
        end
    end;
end;

S1=zeros(size(X));
S2=zeros(size(X));
T1=zeros(size(X));
T2=zeros(size(X));

% Implement Knutsson Eq (29)
for ii=1:NF-1
    S1=S1+q{ii};
    S2=S2+2^((ii+0.5)/2)*q{ii+1};
end;

% Local frequency (Knutsson Eq 29)
r = r0*S2./S1;

% Estimate "certainty" (confidence) (Knutsson Eq 31)
for ii=1:NF-1
    T1=T1+q{ii}.^2;
    T2=T2+ q{ii}.^2.*(2^((ii+0.5)/2)*q{ii+1}./q{ii} -r).^2;
end;
sig2 = T2./T1;
c = 1./(1+sig2);