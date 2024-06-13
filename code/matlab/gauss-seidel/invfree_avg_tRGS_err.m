function [errs,res_errs,ln_errs,res_ln_errs] = invfree_avg_tRGS_err(A,B,X0,T,X_true,X_LN,block_size,weights)
    %record number of row slices
    [~,l,n] = size(A);
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
        %sample column block
        tau_t = randsample(l,block_size);
        A_slice = A(:,tau_t,:);
        I = eye(l,l);
        E = zeros(l,block_size,n);
        E(:,:,1)=I(:,tau_t);


        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_norm = norm(A_slice_t(:),'fro');
        resid = B - tprod(A,X);

        %GS step
        X = X + (weights(t)/A_prod_norm^2)*tprod(tprod(E,A_slice_t),resid);
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
