function [X,its] = QTRK_new(A,B, X0, T, q_value) 
    X = X0; %initialize iterate
    its = {X}; % storing all approximations 
    %iterate
    for t = 1:T
        E = abs(tprod(A, its{t})-B); %error tensor
        Q = quantile(E, q_value, "all");
        corrupted_rows = [];
         % Identify corrupted rows
        for i = 1:size(B, 1)
            if any(E(i, :, :) >= Q, 'all')
                corrupted_rows(end + 1) = i;
            end
        end

        total_rows = 1:size(B,1);
        uncorrupted_rows = setdiff(total_rows, corrupted_rows);

        % Check if there are uncorrupted rows left
        if isempty(uncorrupted_rows)
            %warning('All rows are corrupted at iteration %d. Stopping early.', t);
            its{end + 1} = X; % Store the current X without updating
            continue; % Skip to the next iteration
        end

        %sample row slice from uncorrupted ones
        i_t = randsample(uncorrupted_rows, 1, true);
        A_slice = A(i_t,:,:);
        B_slice = B(i_t,:,:);

        %calculate necessary transformations of A_slice
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice,A_slice_t));
        resid = tprod(A_slice,X) - B_slice;
        temp_subtract = tprod(tprod(A_slice_t,A_prod_inv),resid);

        %display(['Iteration: ', num2str(t), 'Row: ', num2str(i_t), 'Bad rows: ', num2str(corrupted_rows)]);
        X = X - temp_subtract;
        its{end+1} = X;
    end
end