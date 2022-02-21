function [nmse,psnr] = result_evaluation_phase(result,reference,ifabs)
% ------------------------------------------------------
% Usage:
% used to evaluate the quality of reconstruction result
% Input: 
% (1) reconstructed image
% (2) reference image
% ps: the size should be the same
% 
% Output:
% rmse, psnr
% ------------------------------------------------------
% Runke Wang
% 2020/10/06
% ------------------------------------------------------
[x1,y1,NFrame1] = size(result);
[x2,y2,NFrame2] = size(reference);
if NFrame1 ~= NFrame2
    disp('Please input reconstructed results and reference correctly');
    return;
end


for phase = 1:NFrame1
    X = result(:,:,phase);
    ref = reference(:,:,phase);
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
%     psnr(1,phase) = 10*log10(255^2/mse(1,phase));
    psnr(1,phase) = 10*log10(1^2/mse(1,phase));
end

end

function imgNorm = normabs (imgIn)
imgIn = abs(imgIn);
imgNorm = (imgIn - min(imgIn(:))) / (max(imgIn(:))-min(imgIn(:)));
end