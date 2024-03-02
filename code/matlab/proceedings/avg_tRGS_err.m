function [errs,res_errs,ln_errs,res_ln_errs] = avg_tRGS_err(A,B,X0,T,X_true,X_LN,block_size,weights)
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
        update = zeros(size(X));
        for j = 1:block_size
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
            update = update + weights(j)*tprod(tprod(E,tprod(A_prod_inv,A_slice_t)),resid);
        end

        %GS step
        X = X + update;
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
