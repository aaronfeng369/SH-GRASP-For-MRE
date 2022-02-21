function [H, h, g] = ideal_filter_v2(M, N, D0, type)
% Ideal filter
% Input Variables:
% [M, N] - image matrix size
% D0 - filter radius [lowpass limit, highpass limit]
% type - filter type: 'l' - lowpass, 
%                     'h' - highpass, 
%                     'b' - bandpass
%
% Output Variables:
% H - frequency domain filter display
% h - 2D FIR spatial filter
% g - 2D frequency response
% 
% Record of Revision
% Dec.7, 2010=======Yuan Aaron Feng=======Original Code adapted from PVB

% generate frequency space
[f1,f2] = freqspace([M, N], 'meshgrid'); 
% filter radius
r = sqrt(f1.^2 + f2.^2);
% Initializing
H = ones(size(f1));

if type == 'l'
    % low-pass filter
    H(r>D0(1)) = 0;
    figure, imagesc(H); axis image; title('Lowpass filter')
elseif type == 'h'
    % high-pass filter
    H(r<D0(2)) = 0;
    figure, imagesc(H); axis image; title('Highpass filter')
else
    % bandpass filter
    H((r<D0(1)) | (r>D0(2))) = 0;
    figure, imagesc(H); axis image; title('Bandpass filter')
end

% generate filter matrix
h = fsamp2(H);
% frequency response
g = freqz2(h);
% plot frequency response
figure, surf(abs(g)); title('Filter frequency response')

end