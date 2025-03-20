function [X,its] = mtRK_new(A,B,X0,T, q)
% goal:  solve a linear system AX=B, using a masked quantile tensor
% randomized Kaczmarz
% inputs A,B, X0 (initial value for X)
% T = max number of iterations 
% q = quantile level

    X = X0; %initialize iterate
    its = {X}; % storing all approximations 
    rs = randsample(size(A,1),T,true); %random sample of row slices of A, with rep.
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
        cq = quantile(E, q, "all"); %find the quantile value corresponding to q, in E (error tensor) 
        bad_j = []; % indices that are potential corruptions
        for j = 1:size(E,2)
            for k = 1:size(E,3)
                if E(i_t, j, k) > cq && (any(bad_j == j)) == 0 
                    %check for errors that are larger than cq or are
                    %already in bad_j
                   bad_j(end+1) = j;
                end
            end       
        end 
        temp_subtract(:,bad_j,:) = zeros(size(X0, 1),length(bad_j),size(X0, 3)); 
        X = X - temp_subtract;
        its{end+1} = X;
    end
end
