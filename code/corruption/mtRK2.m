function [X,its] = mtRK2(A,B,X0,T, c_sub1, c_sub2)
    %record number of row slices
    m = size(A,1);

    %initialize iterate
    X = X0;
  

    its = {X};
    %rng(3); 
    rs = randsample(m,T,true);
    %iterate
    for t = 1:T
        %sample row slice
       %edit
        %i_t = randsample(m,1);
        i_t = rs(t);
        A_slice = A(i_t,:,:);
        B_slice = B(i_t,:,:);

        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice,A_slice_t));
        resid = tprod(A_slice,X) - B_slice;

        if  any(c_sub1 == i_t) % check if the current row is where the corrupted row 
            p=find(c_sub1 == i_t);
        %if (i_t, c_sub1) == 1
            % then we do partial update
            % TODO
            %tempX = X(:,  [1:c_sub2-1, c_sub2+1:7] ,:) ; 

            temp_subtract = tprod(tprod(A_slice_t,A_prod_inv),resid);
            % assign zeroes to entries that will 
            temp_subtract(:,c_sub2(p),:) = zeros( size(X0, 1),length(c_sub2(p)),size(X0, 3)); 

            X = X - temp_subtract;
        else 
             %RK step
             X = X - tprod(tprod(A_slice_t,A_prod_inv),resid);

        end

        its{end+1} = X;
      
        
    end
end

