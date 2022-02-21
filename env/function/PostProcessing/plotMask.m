function [maskNaN3D, mask03D] = plotMask(temp_image)
% ----------------------------------
% Runke Wang
% Version 2020/06/15
% ----------------------------------
% This function is used for generating mask for the plot the ROI 
% Input:
% Image to be unwrapped
% Output:
% zero mask
% NaN mask
%-----------------------------------
temp_image=squeeze(temp_image);
figure;
subplot(1,2,1)
imagesc(abs(temp_image));axis('off');axis('image');colormap gray;
%custom mask generation
[mask, xi, yi] = roipoly;
hold on; plot(xi,yi,'r.-');hold off
%show the edge of mask
subplot(1,2,2);imagesc(mask);axis('off');axis('image');colormap gray;
maskNaN3D = ones(size(mask));
maskNaN3D(mask == 0) = NaN;
mask03D = ones(size(mask));
mask03D(mask == 0) = 0;
end