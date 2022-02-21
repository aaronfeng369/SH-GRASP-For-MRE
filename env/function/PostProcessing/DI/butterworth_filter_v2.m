function [H, h, g] = butterworth_filter_v2(M, N, D0, n, type)
% Butterworth filter
% Input Variables:
% [M, N] - image matrix size
% D0 - filter radius [lowpass limit, highpass limit]
% n - filter order
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
% Dec.7, 2010=======Yuan Aaron Feng=======Original Code
% Dec.15,2010=======Yuan Aaron Feng=======Add absolute value for surf 

% generate frequency space
[f1,f2] = freqspace([M, N], 'meshgrid'); 
% filter radius
r = sqrt(f1.^2 + f2.^2);

% Butterworth low-pass filter
Hl = 1./(1 + (r./D0(1)).^(2*n));
% Butterworth high-pass filter
Hh = 1./(1 + (D0(2)./r).^(2*n));
% Butterworth bandpass filter
Hb = (1-Hl).*(1-Hh);

if type == 'l'
    H = Hl;
    figure, imagesc(Hl); axis image; title('Lowpass filter')
elseif type == 'h'
    H = Hh;
    figure, imagesc(Hh); axis image; title('Highpass filter')
else
    H = Hb;
    figure, imagesc(Hb); axis image; title('Bandpass filter')
end

% generate filter matrix
h = fsamp2(H);    
% frequency response
g = freqz2(h);  
% plot frequency response
figure, surf(abs(g)); title('Filter frequency response')

%====================
% [U, V] = dftuv(M, N);
% 
% D = sqrt(U.^2 + V.^2);
% 
% % low-pass
% hl = 1./(1 + (D./D0(2)).^(2*n));
% [Hl, f1l, f2l] = freqz2(hl);
% figure, mesh(abs(Hl));
% 
% % high-pass
% hh = 1./(1 + (D0(1)./D).^(2*n));
% [Hh, f1h, f2h] = freqz2(hh);
% % hl2 = 1./(1 + (D./D0(2)).^(2*n));
% % [Hl2, f12, f22] = freqz2(hl2);
% % figure, mesh(abs(Hl2));
% % Hh = 1 - abs(Hl2);
% % figure, mesh(abs(Hh));
% 
% [Hb,f1b,f2b] = freqz2(hh.*hl);
% figure, mesh(abs(Hb))
% 
% 
% figure, mesh(abs(Hl.*Hh))

end