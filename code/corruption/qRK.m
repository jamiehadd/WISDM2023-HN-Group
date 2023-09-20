function [X,its] = qRK(A,B,X0,T)
    %record number of row slices
    m = size(A,1);
    %initialize iterate
    X = X0;
    its = {X};
    %iterate
    for t = 1:T
        %sample row slice
        i_t = randsample(m,1);
        A_slice = A(i_t,:,:);
        B_slice = B(i_t,:,:);
        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice,A_slice_t));
        resid = tprod(A_slice,X) - B_slice;
        resid_mat = reshape(resid, size(resid,2), size(resid,3));
        resid_fro = norm(resid_mat, 'fro');
        %display(['Chosen Fro norm is: ', num2str(resid_fro)]);
        %Now compare against all other residuals
        q = 0;
        for i=1:m 
            if i~= i_t
                A_slice_tmp = A(i,:,:);
                B_slice_tmp = B(i,:,:);
                resid_tmp = tprod(A_slice_tmp,X) - B_slice_tmp;
                resid_tmp_mat = reshape(resid_tmp,size(resid_tmp,2),size(resid_tmp,3));
                resid_tmp_fro = norm(resid_tmp_mat,'fro');
                if resid_fro > resid_tmp_fro
                    q = q+1;
                end
            end
        end
        q = q/m; %Percentile of our selected slice
        if q > 1
            display(['Iteration: ', num2str(t), ' -- Skipping slice ', num2str(i_t), ', q = ', num2str(q)]);
        else
            %RK step
            X = X - tprod(tprod(A_slice_t,A_prod_inv),resid); 
        end
        its{end+1} = X;
    end
end
