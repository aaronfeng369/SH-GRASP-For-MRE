function [result_display,result_display_cat] = post_processing(result,result_mag,K,name,color)
% ------------------------------------------------------
% Usage:
% used to post-procss and display the result of reconstructed image
% Input: 
% (1) result of reconstructed image
% (2) parameters of the motion encoding gradient 
% Output:
% (1) result of post-precessing
% (2) result after concatenation
% ------------------------------------------------------
% Runke Wang
% 2020/10/20
% ------------------------------------------------------
% K=2*pi*(42.58*10^6)*(amp*10^(-3))*((ramptime+flattoptime)*10^(-6));
NFrame = 4;

x_target = 128;
y_target = 128;
scale = [-1.5*10^(-5),1.5*10^(-5)]; % scale for displacement
%scale = [-1.5,1.5]; % scale for phase

%% remove temporal and spatial DC component
result_display_0 = result-meanphase_temporal(result);
for phase = 1:NFrame 
    %crop and rotate
    result_display(:,:,phase) = rot90(crop(result_display_0(:,:,phase),x_target,y_target),-1);
    %remove direct component
    result_display(:,:,phase) = ( (result_display(:,:,phase)-meanphase_spatial(result_display(:,:,phase)))./K) ./2;
end

%% Gausian filter
result_display = DENSE_filter(result_display);
% Change nan to 0
result_display(isnan(result_display))=0;

%% Concate and display
for phase = 1:NFrame 
    result_display_cat(:,(phase-1)*x_target+1:phase*x_target)=result_display(:,:,phase);
end
figure;imshow(result_display_cat,[],'border','tight','initialmagnification','fit');caxis(scale);
colorbar;colormap(color);
set (gcf,'Position',[0,0,x_target*4+100,y_target]);
h = colorbar;h.Label.String = 'Displacement (μm)';

%% image store
for phase = 1:NFrame 
    figure(20+phase);imshow(result_display(:,:,phase),'border','tight','initialmagnification','fit');set(gcf,'Position',[0,0,x_target,y_target]);
    caxis(scale);
    colormap(color)
    saveas((20+phase), [name,num2str(phase),'.bmp']);
end

for i=1:NFrame   
    I=result_display(:,:,i); 
    figure(50);
    imshow(I,[],'border','tight','initialmagnification','fit');caxis(scale);
    colorbar;colormap(color);
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if i==1
        imwrite(I,map,[name,'_phase.gif'],'gif','WriteMode','overwrite', 'Loopcount',inf,'DelayTime',0.5);
    elseif i==NFrame
        imwrite(I,map,[name,'_phase.gif'],'gif','WriteMode','append','DelayTime',0.5);
    else
        imwrite(I,map,[name,'_phase.gif'],'gif','WriteMode','append','DelayTime',0.5);
    end
    close 50
end


end