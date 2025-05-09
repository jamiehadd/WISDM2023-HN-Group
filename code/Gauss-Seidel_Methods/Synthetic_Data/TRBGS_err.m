function [rel_errs,res_errs,ln_errs,res_ln_errs,abs_errs] = TBRGS_err(A,B,X0,T,X_true,X_ln,mu)
%--------------------------------------------------------------------------
% TBRGS_err: Tensor Randomized Block Gauss Seidel ERRor
% [rel_errs,res_errs,ln_errs,res_ln_errs,abs_errs] = TBRGS_err(A,B,X0,T,X_true,X_ln,mu)
% 
% This routine computes different types of errors for the regression 
% problem
%                   X_star = min_X ||AX - B|| 
% using the TRBGS (tensor block randomized Gauss Seidel) algorithm. 
%
% Inputs:
% A         -- measurement matrix
% B         -- constant tensor
% T         -- number of iterations
% X_0       -- initialization of X
% X_true    -- exact solution
% X_ln      -- least norm solution
% mu        -- block size
%
% Outputs
% rel_errs  -- relative error
% res_errs  -- residual error
% ln_errs   -- least norm error
% res_ln_errs- residual least norm error
% abs_errs  -- absolute error
%--------------------------------------------------------------------------
% Authors: WiSDM 2023 Tensor Group (A. Castillo, J. Haddock, Iryna Hartsock, 
%          P. Hoyos, L. Kassab, A. Kryshchenko, K. Larripa, D. Needell, S.
%          Suryanarayanan, K. Yacoubou Djima)
% 
% Last Modified: May 6, 2025
%--------------------------------------------------------------------------

    % Record number of row slices
    [~,l,n] = size(A);
    abs_errs = [];
    rel_errs = [];
    ln_errs = [];
    res_errs = [];
    res_ln_errs = [];

    % Initialize iterate
    X = zeros(size(X0));
    est = X - X_true;
    ln_est = X - X_ln;
    res_est = tprod(A,est);
    res_ln_est = tprod(A,ln_est);
    abs_errs = [abs_errs,norm(est(:))];
    rel_errs = [rel_errs,norm(est(:))/norm(X_true(:))];
    ln_errs = [ln_errs,norm(ln_est(:))/norm(X_ln(:))];
    res_errs = [res_errs,norm(res_est(:))];
    res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];

    % Iterates computation
    for t = 1:T
        % Sample column slice
        tau = randsample(l,mu);
        A_tau = A(:,tau,:);
        I = eye(l,l);
        E = zeros(l,mu,n);
        E(:,:,1)=I(:,tau);

        % Calculate necessary transformations of A_slice
        A_tau_pinv = tpinv(A_tau);
        resid = B - tprod(A,X);

        % BRGS step
        X = X + tprod(E,tprod(A_tau_pinv,resid));
        est = X - X_true;
        ln_est = X - X_ln;
        res_est = tprod(A,est);
        res_ln_est = tprod(A,ln_est);
        
        abs_errs = [abs_errs,norm(est(:))/norm(X_true(:))];
        rel_errs = [abs_errs,norm(est(:))/norm(X_true(:))];
        ln_errs = [ln_errs,norm(ln_est(:))/norm(X_ln(:))];
        res_errs = [res_errs,norm(res_est(:))];
        res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
    end
end
