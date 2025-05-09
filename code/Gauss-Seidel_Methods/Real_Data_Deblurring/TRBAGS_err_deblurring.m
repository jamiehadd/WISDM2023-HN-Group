function [X,res_errs] = TRBAGS_err_deblurring(A,B,X0,T,X_true,block_size,weights)
%--------------------------------------------------------------------------
% TRBAGS_err_deblurring: Average (pseudo-inverse free) Tensor Randomized Block
%             Gauss Seidel ERRor Deblurring
% [X,res_errs] = TRBAGS_err_deblurring(A,B,X_0,T,X_true,mu,w)
%
% This routine performs deblurring on a data matrix using TRBAGS
%
% Inputs:
% A         -- measurement matrix
% B         -- constant tensor
% T         -- number of iterations
% X_0       -- initialization of X
% X_true    -- exact solution
% mu        -- block size
% w         -- weights
%
% Outputs
% X         -- least norm solution
% res_errs  -- residual 
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
    X = X0;
    est = X - X_true;
    res_est = tprod(A,est);
    res_errs = [res_errs,norm(res_est(:))];

    %iterate
    for t = 1:T
        % Sample column block
        tau_t = randsample(l,block_size);
        A_slice = A(:,tau_t,:);
        I = eye(l,l);
        E = zeros(l,block_size,n);
        E(:,:,1)=I(:,tau_t);


        % Salculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_norm = norm(A_slice_t(:),'fro');
        resid = B - tprod(A,X);

        % GS step
        X = X + (weights(t)/A_prod_norm^2)*tprod(tprod(E,A_slice_t),resid);
        est = X - X_true;
        res_est = tprod(A,est);
        res_errs = [res_errs,norm(res_est(:))];
    end
end
