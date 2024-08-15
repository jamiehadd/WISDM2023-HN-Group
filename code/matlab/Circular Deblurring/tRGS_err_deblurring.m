function [X,res_errs] = tRGS_err_deblurring(A,B,X0,T,X_true)
    %record number of row slices
    [~,l,n] = size(A);
    res_errs = [];

    %initialize iterate
    X = X0;
    est = X - X_true;
    res_est = tprod(A,est);

    res_errs = [res_errs,norm(res_est,"fro")/norm(B,"fro")];

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
        res_est = tprod(A,X) - B;
        res_errs = [res_errs,norm(res_est,"fro")/norm(B,"fro")];
    end
end
