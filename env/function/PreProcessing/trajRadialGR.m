function traj = trajRadialGR(Nkx, Nspokes, overSample, dx,dy)
% Function traj_gr is to generate trajactory lines in golden-ratio radial sampling
%
% Input Variables:
% Nkx - read-out sampling points
% Nspokes - sampling radial spokes
% overSample - over sampling ratio in the read-out direction
% dx, dy - kx and ky sampling skipping ratio
%
% Output Variables:
% traj - trajactory
%
% Record of Revisions:
% May-15-2020===Yuan Feng===Re-write for radial sampling
% May-16-2020===Yuan Feng===Add over-sampling ratio

GR = (sqrt(5)-1)/2;  % Golden Ratio
% GR = 2*pi/Nkx;  % Sequenctial
Nkx = overSample*Nkx;
traj = zeros(Nkx,Nspokes);
% t = zeros(kx*2,spokes);  % double oversampling - this is for UIH scanner
theta = (0:Nspokes-1) * GR * pi;

kk = (-Nkx/2:Nkx/2-1)/Nkx;
traj = kk'*[cos(theta)+1i*sin(theta)];

if nargin > 3
    % This part needs to be updated.
%     dx = 0;  % kx delay
%     dy = 0;  % ky delay
    Coef_x = zeros(1,Nspokes);
    Coef_y = zeros(1,Nspokes);
    
    for jj = 1 : Nspokes
        Coef_x(jj) = cos(theta(jj));
        Coef_y(jj) = sin(theta(jj));
        for ii = 1:Nkx % kx*2 - for double oversampling
            kk = -Nkx/2 + (ii - 1); % - kx/2 + (ii - 1)/2 - for double oversampling
            traj(ii,jj) = (kk-dx) * Coef_x(jj) + 1i * (kk-dy) * Coef_y(jj);
        end
    end
    traj = traj/Nkx;
end



