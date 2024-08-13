function [X,its] = mtRK(A,B,X0,T, c_sub1, c_sub2)
    X = X0;  %initialize iterate
    its = {X};
    rs = randsample(size(A,1),T,true);
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
        
        if  any(c_sub1 == i_t) % check if the current row is where the corruption is 
            p=find(c_sub1 == i_t);
            % assign zeroes to entries in X that are affected by the corruptions 
            temp_subtract(:,c_sub2(p),:) = zeros(size(X0, 1),length(c_sub2(p)),size(X0, 3)); 
            X = X - temp_subtract;
        else 
             X = X - tprod(tprod(A_slice_t,A_prod_inv),resid);
        end
        its{end+1} = X; 
    end
end

