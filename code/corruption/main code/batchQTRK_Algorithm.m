function [X,its] = batchQTRK_Algorithm(A,B,X0,K,q,T) 

% Goal: Solve a tensor linear system AX=B using QTRK
% Inputs: A,B,X0 (initial value for unknown tensor X)
% K = max number of iterations 
% q = quantile level
% T = batch size

X = X0; % initialize algorithm 
its = {X}; % store all approximations 
total_rows = 1:size(B,1); % array of row indices

for t = 1:K
    
    % Compute residual tensor
    E = abs(tprod(A, its{t})-B); 
    Q = quantile(E, q, "all");
    
    % Identify uncorrupted rows
    corrupted_rows = [];
    for i = 1:size(B, 1)
        if any(E(i, :, :) > Q, 'all')
            corrupted_rows(end + 1) = i;
        end
    end
    uncorrupted_rows = setdiff(total_rows, corrupted_rows);

    % Check for uncorrupted rows
    if isempty(uncorrupted_rows)
        %disp("Warning: Stopped early at iteration: " + t);
        for h = length(its)+1:K+1
            its{h} = X;
        end
        return; 
    end
  
    resid = tprod(A,X) - B;
    % Sample row slices from uncorrupted ones
    S_t = randsample(uncorrupted_rows, T);
        
    % iterate through sampled rows
    for j = 1:T
        i_t = S_t(j); % current sampled row index
        
        A_slice = A(i_t,:,:);
        B_slice = B(i_t,:,:);
        resid_slice = resid(i_t,:,:);
    
        % Compute the projection
        A_slice_t = tran(A_slice);
        A_prod_inv = tinv(tprod(A_slice,A_slice_t));
        t_proj = tprod(tprod(A_slice_t,A_prod_inv),resid_slice);
        
        % Update Approximation
        X = X - t_proj;
    end

    its{end+1} = X;

end

end