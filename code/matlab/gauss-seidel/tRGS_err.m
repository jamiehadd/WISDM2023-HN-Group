function [X,errs,res_errs,ln_errs,res_ln_errs,norms] = tRGS_err(A,B,X0,T,X_true,X_LN)
    %record number of row slices
    [m,l,n] = size(A);
    errs = [];
    ln_errs = [];
    res_errs = [];
    res_ln_errs = [];
    norms = [];

    %initialize iterate
    X = X0;
    est = X - X_true;
    ln_est = X - X_LN;
    res_est = tprod(A,est);
    res_ln_est = tprod(A,ln_est);
    errs = [errs,norm(est(:))/norm(X_true(:))];
    ln_errs = [ln_errs,norm(ln_est(:))/norm(X_LN(:))];
    res_errs = [res_errs,norm(res_est(:))];
    res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
    norms = [norms,norm(X(:))];

    %iterate
    for t = 1:T
        %sample column slice
        i_t = randsample(l,1);
        A_slice = A(:,i_t,:);
        I = eye(l,l);
        E = zeros(l,1,n);
        E(:,:,1)=I(:,i_t);


        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice_t,A_slice));
        resid = B - tprod(A,X);

        %RGS step
        X = X + tprod(tprod(E,tprod(A_prod_inv,A_slice_t)),resid);
        est = X - X_true;
        ln_est = X - X_LN;
        res_est = tprod(A,X) - B;
        res_ln_est = tprod(A,ln_est);
        errs = [errs,norm(est(:))/norm(X_true(:))];
        ln_errs = [ln_errs,norm(ln_est(:))/norm(X_LN(:))];
        res_errs = [res_errs,norm(res_est(:))];
        res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
        norms = [norms,norm(X(:))];
    end
end
