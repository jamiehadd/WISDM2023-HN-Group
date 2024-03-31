function [X,errs,res_errs,ln_errs,res_ln_errs,norms] = qTRK_err(A,B,X0,T, qv, X_true,X_LN)
    %record number of row slices
    m = size(A,1);
    errs = [];
    ln_errs = [];
    res_errs = [];
    res_ln_errs = [];
    norms = [];

    %initialize iterate
    X = X0;
    est = X - X_true;
    ln_est = X - X_LN;
    res_est = tprod(A,X)-B;
    res_ln_est = tprod(A,ln_est);
    errs = [errs,norm(est(:))/norm(X_true(:))];
    ln_errs = [ln_errs,norm(ln_est(:))/norm(X_LN(:))];
    res_errs = [res_errs,norm(res_est(:))];
    res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
    norms = [norms,norm(X(:))];

    %iterate
    tic
    for t = 1:T
        if mod(t,1000) == 0
            toc
        end
        %sample row slice
        i_t = randsample(m,1);
        A_slice = A(i_t,:,:);
        B_slice = B(i_t,:,:);
        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice,A_slice_t));
        resid = tprod(A_slice,X) - B_slice;
        resid_mat = reshape(resid, size(resid,2), size(resid,3));
        resid_fro = norm(resid_mat, 'fro');
        %display(['Chosen Fro norm is: ', num2str(resid_fro)]);
        %Now compare against all other residuals
        q = 0;
        for i=1:m 
            if i~= i_t
                A_slice_tmp = A(i,:,:);
                B_slice_tmp = B(i,:,:);
                resid_tmp = tprod(A_slice_tmp,X) - B_slice_tmp;
                resid_tmp_mat = reshape(resid_tmp,size(resid_tmp,2),size(resid_tmp,3));
                resid_tmp_fro = norm(resid_tmp_mat,'fro');
                if resid_fro > resid_tmp_fro
                    q = q+1;
                end
            end
        end
        q = q/m; %Percentile of our selected slice
        if q <= qv
            %RK step
            X = X - tprod(tprod(A_slice_t,A_prod_inv),resid); 
        end
        est = X - X_true;
        ln_est = X - X_LN;
        res_est = tprod(A,X)-B;
        res_ln_est = tprod(A,ln_est);
        errs = [errs,norm(est(:))/norm(X_true(:))];
        ln_errs = [ln_errs,norm(ln_est(:))/norm(X_LN(:))];
        res_errs = [res_errs,norm(res_est(:))];
        res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
        norms = [norms,norm(X(:))];
    end
end
