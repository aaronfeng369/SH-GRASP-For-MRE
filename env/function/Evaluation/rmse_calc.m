function [nmse,psnr] = rmse_calc(result,reference,mode)
% ------------------------------------------------------
% Usage:
% used to evaluate the quality of reconstruction result
%
% Input: 
% (1) reconstructed image
% (2) reference image
% (3) mode: 0 - phase images
%           1 - magnitude images
% ps: the size should be the same
% 
% Output:
% rmse, psnr
% ------------------------------------------------------
% Runke Wang
% 2020/10/06
% ------------------------------------------------------
[x1,y1,NFrame1] = size(result);
[~,~,NFrame2] = size(reference);
if NFrame1 ~= NFrame2
    disp('Please input reconstructed results and reference correctly');
    return;
end

for phase = 1:NFrame1
    if mode==0
        X = result(:,:,phase);
        ref = reference(:,:,phase);
        stiffness_max = max(max(X(:)),max(ref(:)));
        X = X/stiffness_max;
        ref = ref/stiffness_max;
    elseif mode ==1
        X = normabs(result(:,:,phase));
        ref = normabs(reference(:,:,phase));
    end
    mse(1,phase) = 0;
    nmse_d(1,phase) = 0;
    num = 0;
    for i=1:x1
        for j=1:y1
            if ~isnan(ref(i,j))
                mse(1,phase) = mse(1,phase)+(ref(i,j)-X(i,j)).^2;
                nmse_d(1,phase) = nmse_d(1,phase)+ref(i,j).^2;
                num=num+1;
            end
        end
    end
    mse(1,phase) = mse(1,phase)/num;
    nmse(1,phase) = sqrt(mse(1,phase)/nmse_d(1,phase))*100;
    rmse(1,phase) = sqrt(mse(1,phase));
    psnr(1,phase) = 10*log10(1^2/mse(1,phase));
end

end

function imgNorm = normabs (imgIn)
imgIn = abs(imgIn);
imgNorm = (imgIn - min(imgIn(:))) / (max(imgIn(:))-min(imgIn(:)));
end