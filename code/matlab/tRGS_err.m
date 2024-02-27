function errs = tRGS_err(A,B,X0,T,X_true)
    %record number of row slices
    [m,l,n] = size(A);
    errs = [];

    %initialize iterate
    X = X0;
    est = X - X_true;
    errs = [errs,norm(est(:))];
    

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
        errs = [errs,norm(est(:))];
    end
end
