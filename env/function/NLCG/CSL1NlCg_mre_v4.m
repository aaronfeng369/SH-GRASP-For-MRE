function [x,cost] = CSL1NlCg_mre_v3(x0,param)
% 
% res = CSL1NlCg(param)
%
% Compressed sensing reconstruction of undersampled k-space MRI data
%
% L1-norm minimization using non linear conjugate gradient iterations
% 
% Given the acquisition model y = E*x, and the sparsifying transform W, 
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda * ||W*x||_1 
%
% Based on the paper: Sparse MRI: The application of compressed sensing for rapid MR imaging. 
% Lustig M, Donoho D, Pauly JM. Magn Reson Med. 2007 Dec;58(6):1182-95.
%
% Ricardo Otazo, NYU 2008
%
%---------------------------------------------------------------------------------------------
% 
% Modified for MRE data
% (1) add temporal fft for regularization term in function: objective and grad
% (2) add two function for fft and ifft in nlcg iteration
% (3) add modified TV transform: mreTV
% (4) add spatial TV
%
% the program finds the x that minimizes the following objective function:
% (1) f(x) = ||E*x - y||^2 + lambda_fft * ||fft(x)||_1
% (2) f(x) = ||E*x - y||^2 + lambda_mreTV * ||mreTV(x)||_1
% (3) f(x) = ||E*x - y||^2 + lambda_mreTV * ||mreTV(x)||_1 + lambda_fft * ||fft(x)||_1 + lambda_sTV * ||sTV(x)||_1
%
% Input: 
% x0 - Initial images
% param -  parameter setting for cs
%
% Output: reconstructed images

% modified by Runke Wang
% 2021/03/10
%---------------------------------------------------------------------------------------------
    fprintf('\n Non-linear conjugate gradient algorithm')
    fprintf('\n ---------------------------------------------\n')

    % Select the mode: half period shift-0 / half period shift-1
    mode = 0;
    
    % starting point
    x=x0;

    % line search parameters
    maxlsiter = 150 ;
    gradToll = 1e-3 ;
    param.l1Smooth = 1e-15;	
    alpha = 0.01;  
    beta = 0.6;
    t0 = 1 ; 
    k = 0;

    % compute g0  = grad(f(x))
    g0 = grad(x,param,mode);
    dx = -g0;

    % iterations
    while(1)
        % backtracking line-search
        f0 = objective(x,dx,0,param,mode);
        t = t0;
        f1 = objective(x,dx,t,param,mode);
        lsiter = 0;
        while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
            lsiter = lsiter + 1;
            t = t * beta;
            f1 = objective(x,dx,t,param,mode);
        end

        if lsiter == maxlsiter
            disp('Error - line search ...');
            return;
        end

        % control the number of line searches by adapting the initial step search
        if lsiter > 2, t0 = t0 * beta;end 
        if lsiter<1, t0 = t0 / beta; end

        % update x
        x = (x + t*dx);

        % print some numbers	
        if param.display
            fprintf(' ite = %d, cost = %f \n',k,f1);
        end

        %conjugate gradient calculation
        g1 = grad(x,param,mode);
        bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
        g0 = g1;
        dx =  - g1 + bk* dx;
        k = k + 1;
        
        cost(k) = f1;

        % stopping criteria (to be improved)
        if (k > param.nite) || (norm(dx(:)) < gradToll), break;end
        
%         figure(200);
%         for phase =1:4
%             subplot(2,2,phase);imshow(angle(x(:,:,phase)),[]);
%             drawnow();
%         end
    end
    return;
end

function res = objective(x,dx,t,param,mode) %**********************************
    %------Basic-----------
    % L2-norm part
    w=param.E*(x+t*dx)-param.y;
    L2Obj=w(:)'*w(:);

    % L1-norm part
    if param.lambda
       w = param.W*(x+t*dx); 
       L1Obj = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
    else
        L1Obj=0;
    end

    % Spatial TV part
    if param.sTVWeight
       w = param.sTV*(x+t*dx); 
       sTVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
    else
        sTVObj=0;
    end
    
    % Wavelet Part
    if param.WaveletWeight
        w = param.wavelet*(x+t*dx);
        WaveletObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
    else
        WaveletObj=0;
    end
    %------End of Basic-----------

    % Temporal fft l1-norm part
    if param.lambda_fft
           w = fft_nlcg(x+t*dx); 
           L1Obj_fft = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
    else
        L1Obj_fft=0;
    end

    % MRE TV l1-norm part
    if param.lambda_mreTV_L1
        w = mreTV(x+t*dx,mode); 
        L1Obj_mreTV = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
    else
        L1Obj_mreTV = 0;
    end

    % objective function
    res=L2Obj + param.lambda*L1Obj + param.lambda_fft*L1Obj_fft + param.lambda_mreTV_L1*L1Obj_mreTV + param.sTVWeight*sTVObj + param.WaveletWeight*WaveletObj;
end

function g = grad(x,param,mode)%***********************************************
    %------Basic-----------
    % L2-norm part
    L2Grad = 2.*(param.E'*(param.E*x-param.y));
    L2Grad(isnan(L2Grad(:)))=0;

    % L1-norm part
    if param.lambda 
        w = param.W*x;
        L1Grad = param.W'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    else
        L1Grad=0;
    end

    % Spatial TV part
    if param.sTVWeight
        w = param.sTV*x;
        sTVGrad = param.sTV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    else
        sTVGrad=0;
    end
    
    % Wavelet Part
    if param.WaveletWeight
        w = param.wavelet*x;
        WaveletGrad = param.wavelet'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    else
        WaveletGrad=0;
    end
    
    %------End of Basic-----------

    % Temporal fft l1-norm part
    if param.lambda_fft
        w = fft_nlcg(x);
        L1Grad_fft = fft_nlcg_adj(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    else
        L1Grad_fft=0;
    end

    % MRE TV l1-norm part
    if param.lambda_mreTV_L1
        w = mreTV(x,mode);
        w_inv = (w.*conj(w)+param.l1Smooth).^(-0.5);
        L1Grad_mreTV = ( x + mreshift(conj(x),mode) + conj(mreshift_adj(mreshift(conj(x),mode),mode) - mreshift(x,mode)) ).*w_inv;
    else
        L1Grad_mreTV = 0;
    end

    % composite gradient
    g=L2Grad + param.lambda*L1Grad + param.lambda_fft*L1Grad_fft + param.lambda_mreTV_L1*L1Grad_mreTV + param.sTVWeight*sTVGrad + param.WaveletWeight*WaveletGrad;
end

% Temporal fft for MRE regularization term
function w = fft_nlcg(x) 
    w = fftshift(fft(x,[],3),3);
end

function w = fft_nlcg_adj(x) 
    w = ifft(ifftshift(x,3),[],3);
end

% New temporal variation transform for MRE regularization term
function w = mreTV(x,mode)
    if mode == 0
        % Half period shift
        y(:,:,3:8) = x(:,:,1:6);
        y(:,:,1:2) = x(:,:,7:8);
        w = x-conj(y);
    elseif mode == 1
        % Whole period shift
        y(:,:,5:8) = x(:,:,1:4);
        y(:,:,1:4) = x(:,:,5:8);
        w = x-y;
    end
end

function y=mreshift(x,mode)
    if mode ==0
        y(:,:,3:8) = x(:,:,1:6);
        y(:,:,1:2) = x(:,:,7:8);
    elseif mode==1
        y(:,:,5:8) = x(:,:,1:4);
        y(:,:,1:4) = x(:,:,5:8);
    end
end

function y=mreshift_adj(x,mode)
    if mode==0
        y(:,:,7:8) = x(:,:,1:2);
        y(:,:,1:6) = x(:,:,3:8);
    elseif mode==1
        y(:,:,5:8) = x(:,:,1:4);
        y(:,:,1:4) = x(:,:,5:8);
    end
end