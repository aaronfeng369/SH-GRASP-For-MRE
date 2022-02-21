%±£´æ¶¯Ì¬Í¼Ïñ
%clear;
%clc;
phase_num = 4;

SamplePath_imwrite =  'D:\UIH read rawdata\gre_radial_read_rawdata\result\20200602\postprocess\'
fileExt_final = '.fig';
for i=1:phase_num
    num_string = num2str(i);
    picname = strcat(SamplePath_imwrite,num_string,fileExt_final);
    %image=imread(picname);
    %imshow(image,[]);
    open(picname);
    set(gcf,'outerposition');% Matlab window maximization
    frame=getframe(gcf); 
    %To make gif files, images must be index indexed images
    im=frame2im(frame);
    [I,map]=rgb2ind(im,20);          
    if i==1
        imwrite(I,map,'DENSE.gif','gif','WriteMode','overwrite', 'Loopcount',inf,'DelayTime',0.5);
    elseif i==phase_num
        imwrite(I,map,'DENSE.gif','gif','WriteMode','append','DelayTime',0.5);
    else
        imwrite(I,map,'DENSE.gif','gif','WriteMode','append','DelayTime',0.5);
    end
    close all
end