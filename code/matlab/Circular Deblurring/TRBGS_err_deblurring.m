function [X,res_errs] = TRBGS_err_deblurring(A,B,X0,T,X_true,block_size)
    %record number of row slices
    [~,l,n] = size(A);
    res_errs = [];

    %initialize iterate
    X = zeros(size(X0));
    est = X - X_true;
    res_est = tprod(A,est);
    res_errs = [res_errs,norm(res_est(:))];

    %iterate
    for t = 1:T
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
        res_est = tprod(A,est);
       
        res_errs = [res_errs,norm(res_est(:))];
    end
end


