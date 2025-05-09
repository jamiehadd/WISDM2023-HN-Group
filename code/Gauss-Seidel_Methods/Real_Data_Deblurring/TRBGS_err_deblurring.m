function [X,res_errs] = TRBGS_err_deblurring(A,B,X0,T,X_true,block_size)
%--------------------------------------------------------------------------
% TRBGS_err_deblurring: Tensor Randomized BlockGauss Seidel ERRor 
%             Deblurring
% [X,res_errs] = TRBGS_err_deblurring(A,B,X_0,T,X_true,mu,w)
%
% This routine performs deblurring on a data matrix using TRBGS
%
% Inputs:
% A         -- measurement matrix
% B         -- constant tensor
% T         -- number of iterations
% X_0       -- initialization of X
% X_true    -- exact solution
% mu        -- block size
%
% Outputs
% X         -- least norm solution
% res_errs  -- residual errors
%--------------------------------------------------------------------------
% Authors: WiSDM 2023 Tensor Group (A. Castillo, J. Haddock, Iryna Hartsock, 
%          P. Hoyos, L. Kassab, A. Kryshchenko, K. Larripa, D. Needell, S.
%          Suryanarayanan, K. Yacoubou Djima)
% 
% Last Modified: May 6, 2025
%--------------------------------------------------------------------------


    % Record number of row slices
    [~,l,n] = size(A);
    res_errs = [];

    % Initialize iterate
    X = zeros(size(X0));
    est = X - X_true;
    res_est = tprod(A,est);
    res_errs = [res_errs,norm(res_est(:))];

    % Iterate
    for t = 1:T
        tau = randsample(l,block_size);
        A_tau = A(:,tau,:);
        I = eye(l,l);
        E = zeros(l,block_size,n);
        E(:,:,1)=I(:,tau);


        % Calculate necessary transformations of A_slice
        A_tau_pinv = tpinv(A_tau);
        resid = B - tprod(A,X);

        % RGS step
        X = X + tprod(E,tprod(A_tau_pinv,resid));
        est = X - X_true;
        res_est = tprod(A,est);
        res_errs = [res_errs,norm(res_est(:))];
    end
end


