function mQTRK_QTRK_plots(tdims, num_corrupt_array, k_array, q_array, num_trials, num_its, cor_size)

%% Code to generate QTRK plots %%

% Create a new folder to the save figure
rnumb = num2str(randi([1, 1000]));
folderName = ['Exp_',num2str(rnumb)] ;
if ~exist(folderName, 'dir') % check it does not exist
    mkdir(folderName);
end

disp("Experiment Folder:")
disp(rnumb);

%% Running Algorithm with Different Combination of Parameters

% Number of parameters
param1 = length(num_corrupt_array);
param2 = length(k_array);

% Initialize variables to determine y-axis limits
allY_Q = [];
allY_m = [];

% Loop through each combination of parameters to compute aggregated results
for i = 1:param1

    % Define current parameters
    num_corrupt = num_corrupt_array(i);
    q = q_array(i);

    for j = 1:param2

        % Define current parameters
        k = k_array(j);

        % Run QTRK
        [errs_matrix_Q,errs_matrix_m] = mQTRK_QTRK_trials(tdims, num_corrupt, q, k, cor_size, num_trials, num_its);

        % Aggregate results
        % --- >
        % Example: Calculate median errors across trials
        median_errs_Q = median(errs_matrix_Q);
        median_errs_m = median(errs_matrix_m);

        % Collect y-values for later use
        allY_Q = [allY_Q; median_errs_Q];
        allY_m = [allY_m; median_errs_m];

    end
end

%% Enforcing same y-axis range for all semi-log plots

% Determine global y-axis limits
yMin = min([allY_Q(:); allY_m(:)]);
yMax = max([allY_Q(:); allY_m(:)]);

% Choice of markers, colors, and lines for plotting
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.2, 0.2, 0.2]};
markers = {'o', '*', '+'};
lineStyles = {'-', ':', '-.'};
legendLabels = cell(1, 2*param2);

% Initialize index counter for accessing allY
index = 1;

% Generate subplots of different parameter combinations
for i = 1:param1 %num_corrupt

    individual_fig = figure;

    hold on;
    for j = 1:param2 % k

        % Plotting style
        lineStyle = lineStyles{j};
        marker = markers{j};

        median_errs_Q = allY_Q(index,:);
        median_errs_m = allY_m(index,:);

        % Plotting 
        plot(1:100:num_its+1, median_errs_Q(1:100:num_its+1), 'Color', colors{1}, 'Marker', marker, 'MarkerSize', 12, 'LineStyle', lineStyles{1}, 'LineWidth', 3);
        plot(1:100:num_its+1, median_errs_m(1:100:num_its+1), 'Color', colors{2}, 'Marker', marker, 'MarkerSize', 12, 'LineStyle', lineStyles{2}, 'LineWidth', 3);

        % Update the index for the next entry
        index = index + 1;
    end

    ind_ll = 1;
    m = tdims(4);
    for ll = 1:2*param2
        if mod(ll, 2) == 1
            legendLabels{ll} = sprintf('QTRK %.2f', k_array(ind_ll)/m); %$\\beta_{\\text{row}}$ =
        else
            legendLabels{ll} =  sprintf('mQTRK %.2f', k_array(ind_ll)/m); %$\\beta_{\\text{row}}$ =
            ind_ll = ind_ll +1;
        end
    end

    % Set the y-axis to a logarithmic scale
    set(gca, 'YScale', 'log');

    % Set font size for tick labels and texts
    set(gca, 'FontSize', 22); %tick labels
    xlabel('Iteration', 'Interpreter','latex', FontSize=24);
    ylabel('Relative Error', 'Interpreter','latex', FontSize=24);

    % Set x- and y-axis limits
    ylim([yMin, yMax]);
    xlim([0, num_its]);
    legend(legendLabels,'Interpreter','latex', 'FontSize', 22);

    % Save individual subplots
    figFileName = fullfile(folderName, sprintf('subplot_%d.fig', num_corrupt_array(i)));
    savefig(individual_fig, figFileName);
    
    % Save as PNG
    set(gcf, 'Position', [100, 100, 600, 400]);  % [left, bottom, width, height]
    pngFileName = fullfile(folderName, sprintf('subplot_%d.png', num_corrupt_array(i)));
    print(gcf, pngFileName, '-dpng', '-r300');  % Adjust resolution as needed

    close(individual_fig);
    hold off;

end

close all

end