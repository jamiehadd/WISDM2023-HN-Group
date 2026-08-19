function [X,its,times] = blockQTRK_Algorithm_runtime(A,B,X0,K,q,T) 

% Goal: Solve a tensor linear system AX=B using Block QTRK
% Inputs: A,B,X0 (initial value for unknown tensor X)
% K = max number of iterations 
% q = quantile level
% T = block size

X = X0; % initialize algorithm 
its = {X}; % store all approximations 
times = zeros(1, K+1);
total_rows = 1:size(B,1); % array of row indices

t0 = tic;

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
        elapsed = toc(t0);
        for h = length(its)+1:K+1
            its{h} = X;
            times(h) = elapsed;
        end
        return;
    end
  
    resid = tprod(A,X) - B;

    % Sample row slices from uncorrupted ones
    S_t = randsample(uncorrupted_rows, T);

    A_slice = A(S_t,:,:);
    resid_slice = resid(S_t,:,:);
    
    % Compute the projection
    A_slice_t = tran(A_slice);
    A_prod_inv = tinv(tprod(A_slice,A_slice_t));
    t_proj = tprod(tprod(A_slice_t,A_prod_inv),resid_slice);
    
    % Update Approximation
    X = X - t_proj;
    its{end+1} = X;

    times(length(its)) = toc(t0); % cumulative wall-clock time (in seconds)
end

end