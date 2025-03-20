function errs_matrix = QTRK_trials(tdims, num_corrupt, q, k, cor_size, num_trials, num_its)


%% PARAMETERS
% tdims (list): list of dimentions [l,p,n,m] where
% the unkown tensor has shape(l,p,n),
% the left-hand measurement tensor A has shape (m,l,n),
% the right-hand measurement tensor B has shape (m,p,n)

% num_corrupt (integer): actual number of corruptions

% q (0<q<1): quantile value

% k (integer): number of corrupted rows in B that we restrict to (0 < k <= m)

% cor_size(list): [mean_corrupt, deviation_corrupt]

% num_trials (integer>0): number of trials of the experiments for robustness of results

% num_its (integer>0): number of iterations for the iterative algorithm

%%

% Pull given dimensions
l = tdims(1);
p = tdims(2);
n = tdims(3);
m = tdims(4);

% Pull corruption size distribution
mean_corrupt = cor_size(1);
deviation_corrupt = cor_size(2);

% Initialization
errs_matrix = zeros(num_trials, num_its+1);

% Perform trials
for t = 1:num_trials

    %generate tensors
    A = abs(randn(m,l,n));
    X_true = abs(randn(l,p,n));
    B = tprod(A,X_true);

    % Generate random corruption values for each trial
    corruption_values = mean_corrupt + deviation_corrupt*randn(num_corrupt, 1);

    % Generate k random row indices to corrupt
    corrupt_rows = randsample(m, k,false)';

    all_inds = [];

    for i = corrupt_rows
        for j = 1:p
            for k = 1:n
                all_inds = [all_inds; i,j,k];
            end
        end
    end

    corrupt_entries_inds = randsample(size(all_inds,1),num_corrupt,false);

    % Distribute num_corrupt corruptions uniformly across the k rows
    for i = 1:num_corrupt
%         % Select a random row from the chosen rows
%         row_idx = corrupt_rows(randsample(k, 1),false);
% 
%         % Randomly select indices for the other dimensions
%         col_idx = randsample(p, 1);
%         depth_idx = randsample(n, 1);
        %all_inds(corrupt_entries_inds(i),:)
        inds = all_inds(corrupt_entries_inds(i),:);
        row_idx = inds(1);
        col_idx = inds(2);
        depth_idx = inds(3);

        % Apply the corruption value
        B(row_idx, col_idx, depth_idx) = B(row_idx, col_idx, depth_idx) + corruption_values(i);
    end

    % Run QTRK
    [~, its] = QTRK_Algorithm(A,B, randn(l,p,n), num_its, q); % Adjust q as needed

    % Record errors for the current trial
    for j = 1:num_its + 1
        est = its{j} - X_true; % Adjust X_true according to your setup
        errs_matrix(t,j) = norm(est(:))/norm(X_true(:)); % Relative Frobenius error to true solution
    end
end
end