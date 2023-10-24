function [X,its] = mtRK_new(A,B,X0,T, ncorr)
    X = X0; %initialize iterate
    its = {X};
    rs = randsample(size(A,1),T,true);
    %iterate
    for t = 1:T
        %sample row slice
        i_t = rs(t);
        A_slice = A(i_t,:,:);
        B_slice = B(i_t,:,:);
        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice,A_slice_t));
        resid = tprod(A_slice,X) - B_slice;
        temp_subtract = tprod(tprod(A_slice_t,A_prod_inv),resid);
        E = abs(tprod(A, its{t})-B); %error tensor
        q = 1 - ncorr/(size(E,1)*size(E,2)*size(E,3)); %proportion of non-corrupted entries
        cq = quantile(E, q, "all"); 
        bad_j = [];
        for j = 1:size(E,2)
            for k = 1:size(E,3)
                if E(i_t, j, k) > cq && (any(bad_j == j)) == 0
                   bad_j(end+1) = j;
                end
            end       
        end 
        temp_subtract(:,bad_j,:) = zeros(size(X0, 1),length(bad_j),size(X0, 3)); 
        X = X - temp_subtract;
        its{end+1} = X;
    end
end
