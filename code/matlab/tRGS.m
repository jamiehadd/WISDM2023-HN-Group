function [X,its] = tRGS(A,B,X0,T)
    %record number of row slices
    [m,l,n] = size(A);

    %initialize iterate
    X = X0;
    its = {X};
    

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
        its{end+1} = X;
    end
end
