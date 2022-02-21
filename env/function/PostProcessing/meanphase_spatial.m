function mean_phase = meanphase_spatial(temp_image)
% ----------------------------------
% This function is used for calculate the mean phase for the phase image of MRE 
% Input:
% Image with DC component (phase offset)
% Output:
% mean phase of the ROI 
% 
% Record of Revisions
% Jun-15-2020===Runke Wang===Original Code
% Apr-10-2021===Yuan Feng===Minor modifications
%-----------------------------------
temp_image=squeeze(temp_image);
% [NFe,NPe] = size(temp_image);
mean_phase = nanmean(temp_image(:));
% num=0;
% phase_offset = 0;
% for x=1:NFe
%     for y=1:NPe
%         if ~isnan(temp_image(x,y))
%             phase_offset=phase_offset+temp_image(x,y);
%             num=num+1;
%         end
%     end
% end
% mean_phase=phase_offset/num;
