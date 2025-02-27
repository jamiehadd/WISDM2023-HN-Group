function [X,res_errs] = TRBAGS_err_deblurring(A,B,X0,T,X_true,block_size,weights)
    %record number of row slices
    [~,l,n] = size(A);
    res_errs = [];

    %initialize iterate
    X = X0;
    est = X - X_true;
    res_est = tprod(A,est);
    res_errs = [res_errs,norm(res_est(:))];

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
        res_est = tprod(A,est);
        res_errs = [res_errs,norm(res_est(:))];
    end
end
