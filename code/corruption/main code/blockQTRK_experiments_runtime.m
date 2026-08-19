%% BlockQTRK Experiments %%

% This script:
% Generates synthetic corrupted tensor linear systems,
% Runs block QTRK, and
% Plots relative error vs. iteration and relative error vs. runtime (s)

clear; clc; close all;
addpath(genpath(pwd));
addpath('../tproduct toolbox 2.0 (transform)/')
%warning('off','all')

tic

% Shuffle state
rngState = rng('shuffle');

% Create a new folder to save results
exp_id = randi([1, 9999]);
folderName = ['exp-', num2str(exp_id)];
if ~exist(folderName, 'dir') % check it does not exist
    mkdir(folderName);
end

% Save random state
save(fullfile(folderName, 'rngState.mat'), 'rngState');

% Define parameters to generate syntehtic data (see hyperparameters.txt)

l = 5; p = 4; n = 10; m = 25;
block_size_set = 1:5; % this is to experiment with different block sizes

q = 0.975;
betarow = 0.2;

num_corrupt = round((1-q)*(m*p*n));
num_rowcorrupt = round(betarow*m);
mean_corrupt = 100;
deviation_corrupt = 20;

num_trials = 150;
num_its = 2000;

disp("Experiment Folder: " + exp_id);
disp("Dims: l = " + l + ", p = " + p + ", n = " + n + ", m = " + m);
disp("Block sizes: " + mat2str(block_size_set));
disp("q = " + q);

% Initialization
b_idx = 1;
median_errs_matrix = zeros(length(block_size_set),num_its+1);
median_time_matrix = zeros(length(block_size_set),num_its+1);


for block_size = block_size_set % sweep block sizes
    disp("Running with block size: " + block_size)
    % Initialization
    errs_matrix = zeros(num_trials, num_its+1);
    time_matrix = zeros(num_trials, num_its+1);

    % Perform trials
    for k = 1:num_trials

        %generate tensors
        A = randn(m,l,n);
        X_true = randn(l,p,n);
        B = tprod(A,X_true);
        X0 = randn(l,p,n);

        % Generate random corruption values for each trial
        corruption_values = normrnd(mean_corrupt, deviation_corrupt, [num_corrupt, 1]);

        % Generate k random row indices to corrupt
        corrupt_rows = randsample(m, num_rowcorrupt, false);

        % Distribute num_corrupt corruptions uniformly across the k rows
        for i = 1:num_corrupt
            % Select a random row from the chosen rows
            row_idx = corrupt_rows(randsample(num_rowcorrupt, 1));

            % Randomly select indices for the other dimensions
            col_idx = randsample(p, 1);
            depth_idx = randsample(n, 1);

            % Apply the corruption value
            B(row_idx, col_idx, depth_idx) = B(row_idx, col_idx, depth_idx) + corruption_values(i);
        end

        % Run Algorithm (instrumented copy that also returns per-iteration runtime)
        [~, its, times] = blockQTRK_Algorithm_runtime(A,B, X0, num_its, q, block_size); % Adjust q as needed

        % Record errors and elapsed time for the current trial
        for j = 1:num_its + 1
            est = its{j} - X_true; % Adjust X_true according to your setup
            errs_matrix(k,j) = norm(est(:))/norm(X_true(:)); % Relative Frobenius error to true solution
        end
        time_matrix(k,:) = times;
    end
    median_errs_matrix(b_idx,:) = median(errs_matrix);
    median_time_matrix(b_idx,:) = median(time_matrix);
    b_idx = b_idx + 1;
end
disp("Completed all computations.")


%% Plot styling shared by both figures below

disp("Generating plots.")
% Choice of markers, colors, and lines for plotting
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.2, 0.2, 0.2], [0.466, 0.674, 0.188], [0.4940, 0.1840, 0.5560]};
lineStyles = {'-', '--',':', '-.'};

% Build legend labels directly from block_size_set
num_blocks = length(block_size_set);
legendLabels = cell(1, num_blocks);
for b = 1:num_blocks
    legendLabels{b} = strcat('$T = ', num2str(block_size_set(b)), '$');
end
markers = {'o', 'x', 's', '^', 'd', 'v'};

%% Plot 1: Relative error vs. iteration

individual_fig = figure;
hold on

for b = 1:num_blocks
    median_errs = median_errs_matrix(b,:);
    styleIdx = mod(b-1, length(lineStyles)) + 1;
    if b > length(lineStyles)
        markerIdx = mod(b - length(lineStyles) - 1, length(markers)) + 1;
        plot(1:100:num_its+1, median_errs(1:100:num_its+1), 'Color', colors{b}, 'Marker', markers{markerIdx}, 'MarkerSize', 12, 'LineStyle', lineStyles{styleIdx}, 'LineWidth', 5);
    else
        plot(1:num_its+1, median_errs, 'Color', colors{b}, 'LineStyle', lineStyles{styleIdx}, 'LineWidth', 5);
    end
end

% Set the y-axis to a logarithmic scale
set(gca, 'YScale', 'log');

% Set font size for tick labels and texts
set(gca, 'FontSize', 34); %tick labels
xlabel('Iteration', 'interpreter','latex', 'FontSize', 38);
ylabel('Relative Error', 'interpreter','latex', 'FontSize', 38);

legend(legendLabels, 'Interpreter','latex', 'FontSize', 34);

% Set x-axis limits
xlim([0, num_its]);

set(gcf, 'Position', [100, 100, 900, 600]);  % [left, bottom, width, height]

% Save figure (.fig)
figFileName = fullfile(folderName, ['blockQTRK_exp_', num2str(exp_id), '_vs_iteration.fig']);
savefig(individual_fig, figFileName);

% Save figure (.png)
pngFileName = fullfile(folderName, ['blockQTRK_exp_', num2str(exp_id), '_vs_iteration.png']);
print(gcf, pngFileName, '-dpng', '-r300');  % Adjust resolution as needed

% Save figure (.eps)
epsFileName = fullfile(folderName, ['blockQTRK_exp_', num2str(exp_id), '_vs_iteration.eps']);
print(gcf, epsFileName, '-depsc', '-r300');

close(individual_fig);
hold off
disp("Plot vs. iteration saved.")

%% Plot 2: Relative error vs. runtime

individual_fig = figure;
hold on

for b = 1:num_blocks
    median_errs = median_errs_matrix(b,:);
    median_times = median_time_matrix(b,:);
    styleIdx = mod(b-1, length(lineStyles)) + 1;
    if b > length(lineStyles)
        markerIdx = mod(b - length(lineStyles) - 1, length(markers)) + 1;
        plot(median_times(1:100:end), median_errs(1:100:end), 'Color', colors{b}, 'Marker', markers{markerIdx}, 'MarkerSize', 12, 'LineStyle', lineStyles{styleIdx}, 'LineWidth', 5);
    else
        plot(median_times, median_errs, 'Color', colors{b}, 'LineStyle', lineStyles{styleIdx}, 'LineWidth', 5);
    end
end

% Set the y-axis to a logarithmic scale
set(gca, 'YScale', 'log');

% Set font size for tick labels and texts
set(gca, 'FontSize', 34); %tick labels
xlabel('Runtime (s)', 'interpreter','latex', 'FontSize', 38);
ylabel('Relative Error', 'interpreter','latex', 'FontSize', 38);

legend(legendLabels, 'Interpreter','latex', 'FontSize', 34);

% Set x-axis limits to cover the slowest block size's full runtime
xlim([0, max(median_time_matrix(:))]);

set(gcf, 'Position', [100, 100, 900, 600]);  % [left, bottom, width, height]

% Save figure (.fig)
figFileName = fullfile(folderName, ['blockQTRK_exp_', num2str(exp_id), '_vs_runtime.fig']);
savefig(individual_fig, figFileName);

% Save figure (.png)
pngFileName = fullfile(folderName, ['blockQTRK_exp_', num2str(exp_id), '_vs_runtime.png']);
print(gcf, pngFileName, '-dpng', '-r300');  % Adjust resolution as needed

% Save figure (.eps)
epsFileName = fullfile(folderName, ['blockQTRK_exp_', num2str(exp_id), '_vs_runtime.eps']);
print(gcf, epsFileName, '-depsc', '-r300');

close(individual_fig);
hold off
disp("Plot vs. runtime saved.")

%% Save experiment info

tmr = toc;
dt = datetime;
disp("Wall-clock time (in sec): " + tmr)

filePath = fullfile(folderName, 'parameters.txt');
discp = fopen(filePath, 'w');
fprintf(discp, "Date and Time of Experiment: %s\n", dt);
fprintf(discp, "Experiment Folder: %d\n", exp_id);
fprintf(discp, "Dims: l = %d, p = %d, n = %d, m = %d\n", l, p, n, m);
fprintf(discp, "Block sizes: %s\n", mat2str(block_size_set));
fprintf(discp, "Trials: %d, Iters: %d\n", num_trials, num_its);
fprintf(discp, "q = %.5f\n", q);
fprintf(discp, "num_corrupt = %d, num_rowcorrupt = %d\n", num_corrupt, num_rowcorrupt);
fprintf(discp, "mean_corrupt = %.3f, deviation_corrupt = %.3f\n", mean_corrupt, deviation_corrupt);
fprintf(discp, "Wall-clock time (in sec): %.3f\n", tmr);
fclose(discp);
close all
