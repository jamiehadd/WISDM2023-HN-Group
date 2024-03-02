function [errs,res_errs,ln_errs,res_ln_errs] = block_tRGS_err(A,B,X0,T,X_true,X_LN,block_size)
    %record number of row slices
    [m,l,n] = size(A);
    errs = [];
    ln_errs = [];
    res_errs = [];
    res_ln_errs = [];

    %initialize iterate
    X = X0;
    est = X - X_true;
    ln_est = X - X_LN;
    res_est = tprod(A,est);
    res_ln_est = tprod(A,ln_est);
    errs = [errs,norm(est(:))];
    ln_errs = [ln_errs,norm(ln_est(:))];
    res_errs = [res_errs,norm(res_est(:))];
    res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];

    %iterate
    for t = 1:T
        %sample column slice
        tau = randsample(l,block_size);
        A_tau = A(:,tau,:);
        I = eye(l,l);
        E = zeros(l,block_size,n);
        E(:,:,1)=I(:,tau);


        %calculate necessary transformations of A_slice
        A_tau_pinv = tpinv(A_tau);
        resid = B - tprod(A,X);

        %RGS step
        X = X + tprod(E,tprod(A_tau_pinv,resid));
        est = X - X_true;
        ln_est = X - X_LN;
        res_est = tprod(A,est);
        res_ln_est = tprod(A,ln_est);
        errs = [errs,norm(est(:))];
        ln_errs = [ln_errs,norm(ln_est(:))];
        res_errs = [res_errs,norm(res_est(:))];
        res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
    end
end
