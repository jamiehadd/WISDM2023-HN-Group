function [X,outer_errs,outer_res_ln_errs,inner_errs,inner_res_ln_errs] = FacTRBAGS_err(U,V,Y,X_0,X_true,X_ln,T,outb,inb,w_1,w_2)
%--------------------------------------------------------------------------
% FacTRBAGS_err: FACtorized Tensor Randomized Block Average 
%                Gauss Seidel ERRor
%
% This routine computes different types of errors for the regression 
% problem
%                       X_star = min_X ||AX - B||, 
% where A is a tensor in factored form A = UV. The aglorithm used is the
% FactRBGS (pseudoinverse free). The errors are: 
%
% Inputs:
% U, V  -- factors of the measurement tensor A
% B     -- constant tensor
% X_0   -- initialization of X
% X_true-- exact solution
% X_ln  -- least norm solution
% out_b  -- outer block size
% in_b   -- inner block size
% w_1   -- weights for inner block
% w_2   -- weights for outer block
%
% Outputs
% X     -- least norm solution
% outer/inner_errs: Absolute error on the outer/inner systems
% outer/inner_res_ln_errs: Residual error on the outer/inner systems
%--------------------------------------------------------------------------
% Authors: WiSDM 2023 Tensor Group (A. Castillo, J. Haddock, Iryna Hartsock, 
%          P. Hoyos, L. Kassab, A. Kryshchenko, K. Larripa, D. Needell, S.
%          Suryanarayanan, K. Yacoubou Djima)
% 
% Last Modified: April 28, 2025
%--------------------------------------------------------------------------

 % Record tensors sizes and set errors sizes
    [~,l,n] = size(U);
    [~,l_2,~] = size(V);
    outer_errs = [];
    inner_errs = [];
    outer_res_ln_errs = [];
    inner_res_ln_errs = [];
    norms = [];

    % Initialize iterates
    X = zeros(size(X_0));
    Z = tprod(V,X_0);
    Z_true = tprod(V,X_true);
    est = X - X_true;
    est_Z = Z - Z_true;

    % Initialize errors
    ln_est = X - X_ln;
    inner_res_est = tprod(V,ln_est);
    res_ln_est = tprod(U,tprod(V,ln_est));
    outer_errs = [outer_errs,norm(est(:))/norm(X_true(:))];
    inner_errs = [inner_errs,norm(est_Z(:))/norm(Z_true(:))];
    inner_res_ln_errs = [inner_res_ln_errs, norm(inner_res_est(:))];
    outer_res_ln_errs = [outer_res_ln_errs,norm(res_ln_est(:))];
    norms = [norms,norm(X(:))];

    % Compute iterates and errors
    for t = 1:T

        %-------- Z-Update for inner System --------
        % Sample column slice from U to update Z in terms of Y residuals
        i_t = randsample(l,outb,false);
        U_slice = U(:,i_t,:);
        I = eye(l,l);
        E = zeros(l,outb,n);
        E(:,:,1)=I(:,i_t);

        % Calculate necessary transformations of U_slice
        U_slice_t = tran(U_slice);
        U_prod_norm = norm(U_slice_t(:),'fro');
        resid = Y - tprod(U,Z);

        % Pseudoinverse-free RGS step
        Z = Z + (w_1(t)/U_prod_norm^2)*tprod(tprod(E,U_slice_t),resid);

        %-------- X-Update for outer system --------
        % sample column slice from V to update X in terms of Z resiuals
        i_t = randsample(l_2,inb,false);
        V_slice = V(:,i_t,:);
        I = eye(l_2,l_2);
        E = zeros(l_2,inb,n);
        E(:,:,1)=I(:,i_t);

        % Calculate necessary transformations of U_slice
        V_slice_t = tran(V_slice);
        V_prod_norm = norm(V_slice_t(:),'fro');
        resid = Z - tprod(V,X);

        % Pseudoinverse-free RGS step
        X = X + (w_2(t)/V_prod_norm^2)*tprod(tprod(E,V_slice_t),resid);
        est = X - X_true;
        est_Z = Z - Z_true;
        ln_est = X - X_ln;
        inner_res_est = tprod(V,ln_est);

        res_ln_est = tprod(U,tprod(V,ln_est)); 
        outer_errs = [outer_errs,norm(est(:))/norm(X_true(:))];
        inner_errs = [inner_errs,norm(est_Z(:))/norm(Z_true(:))]; 
        inner_res_ln_errs = [inner_res_ln_errs, norm(inner_res_est(:))]; 
        outer_res_ln_errs = [outer_res_ln_errs,norm(res_ln_est(:))];
        norms = [norms,norm(X(:))];

    end
end
