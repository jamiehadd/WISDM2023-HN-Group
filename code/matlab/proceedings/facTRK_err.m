function [X,errs,res_errs,ln_errs,res_ln_errs,norms] = facTRK_err(U,V,Y,X_0,Z_true,X_true,X_LN,T)
    %record number of row slices
    m = size(U,1);
    m_2 = size(V,1);
    errs = [];
    ln_errs = [];
    res_errs = [];
    res_ln_errs = [];
    norms = [];

    %initialize iterate
    X = X_0;
    Z = tprod(V,X_0);
    est = X - X_true;
    ln_est = X - X_LN;
    res_est = tprod(U,tprod(V,est));
    res_ln_est = tprod(U,tprod(V,ln_est));
    errs = [errs,norm(est(:))/norm(X_true(:))];
    ln_errs = [ln_errs,norm(ln_est(:))/norm(X_LN(:))];
    res_errs = [res_errs,norm(res_est(:))];
    res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
    norms = [norms,norm(X(:))];

    %iterate
    for t = 1:T
        %sample row slice from U to update Z in terms of Y residuals
        i_t = randsample(m,1);
        U_slice = U(i_t,:,:);
        Y_slice = Y(i_t,:,:);

        %calculate necessary transformations of U_slice
        U_slice_t = tran(U_slice);
        U_prod_inv = tinv(tprod(U_slice,U_slice_t));
        resid_y = tprod(U_slice,Z) - Y_slice;

        %RK step 
        Z = Z - tprod(tprod(U_slice_t,U_prod_inv),resid_y);

        %sample row slice from V to update X in terms of Z residuals
        i_t = randsample(m_2,1);
        V_slice = V(i_t,:,:);
        Z_slice = Z(i_t,:,:);


        %calculate necessary transformations of V_slice
        V_slice_t = tran(V_slice);
        V_prod_inv = tinv(tprod(V_slice,V_slice_t));
        resid_z = tprod(V_slice,X) - Z_slice;

        %RK step 
        X = X - tprod(tprod(V_slice_t,V_prod_inv),resid_z);

        est = X - X_true;
        ln_est = X - X_LN;
        res_est = tprod(U,tprod(V,X)) - Y;
        res_ln_est = tprod(U,tprod(V,ln_est));
        errs = [errs,norm(est(:))/norm(X_true(:))];
        ln_errs = [ln_errs,norm(ln_est(:))/norm(X_LN(:))];
        res_errs = [res_errs,norm(res_est(:))];
        res_ln_errs = [res_ln_errs,norm(res_ln_est(:))];
        norms = [norms,norm(X(:))];
    end
end
