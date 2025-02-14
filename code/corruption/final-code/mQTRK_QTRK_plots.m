function mQTRK_QTRK_plots(tdims, num_corrupt_array, k_array, q_array, num_trials, num_its, cor_size, corr_option)

%% Code to generate QTRK plots %%
tic

% Create a new folder to the save figure
rnumb = num2str(randi([1, 1000]));
folderName = ['Exp_',num2str(rnumb)] ;
if ~exist(folderName, 'dir') % check it does not exist
    mkdir(folderName);
end

% Pull given dimensions
l = tdims(1);
p = tdims(2);
n = tdims(3);
m = tdims(4);

disp("Experiment Folder: " + rnumb);
disp("Corr. dist. = " + cor_size);
disp("beta-row = " + k_array/m);
disp("beta = " + num_corrupt_array/(m*p*n));
disp("q = " + q_array);


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

       disp(i + " out of " + param1 + " and " + j + " out of " + param2)

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

% Generate subpfig of different parameter combinations
for i = 1:param1 %num_corrupt

    individual_fig = figure;

    hold on;
    for j = 1:param2 % k

        median_errs_Q = allY_Q(index,:);
        median_errs_m = allY_m(index,:);

        % Median Plotting 
        plot(1:100:num_its+1, median_errs_Q(1:100:num_its+1), 'Color', colors{1}, 'Marker', markers{j}, 'MarkerSize', 12, 'LineStyle', lineStyles{1}, 'LineWidth', 4);
        plot(1:100:num_its+1, median_errs_m(1:100:num_its+1), 'Color', colors{2}, 'Marker', markers{j}, 'MarkerSize', 12, 'LineStyle', lineStyles{2}, 'LineWidth', 4);

        % Update the index for the next entry
        index = index + 1;
    end

    ind_ll = 1;
    for ll = 1:2*param2
        if mod(ll, 2) == 1
            legendLabels{ll} = strcat('QTRK ' , ' $\beta_{\mbox{row}} = $ ', num2str(round(k_array(ind_ll)/m,3)));
        else
            legendLabels{ll} =  strcat('mQTRK ' , ' $\beta_{\mbox{row}} = $ ', num2str(round(k_array(ind_ll)/m,3)));
            ind_ll = ind_ll +1;
        end
    end

    % Set the y-axis to a logarithmic scale
    set(gca, 'YScale', 'log');

    % Set font size for tick labels and texts
    set(gca, 'FontSize', 32); %tick labels
    xlabel('Iteration', 'Interpreter','latex', FontSize=36);
    ylabel('Relative Error', 'Interpreter','latex', FontSize=36);

    % Set x- and y-axis limits
    ylim([yMin, yMax]);
    xlim([0, num_its]);
    
    if  (i == 1) && (corr_option == "large") % show legend only once 
        legend(legendLabels,'Interpreter','latex', 'FontSize', 30);
    end

    % Save individual subfigures
    corr_option = char(corr_option);
    figFileName = fullfile(folderName, ['COMP_', corr_option, '_exp_', num2str(rnumb), '_subfig_', num2str(num_corrupt_array(i)),'.fig']);
    savefig(individual_fig, figFileName);
    
    % Save as PNG
    set(gcf, 'Position', [100, 100, 600, 400]);  % [left, bottom, width, height]
    pngFileName = fullfile(folderName, ['COMP_', corr_option, '_exp_', num2str(rnumb), '_subfig_', num2str(num_corrupt_array(i)),'.png']);
    print(gcf, pngFileName, '-dpng', '-r300');  % Adjust resolution as needed

    close(individual_fig);
    hold off;

end
tmr = toc;
disp("Wall-clock time (in sec): "  + tmr)

filePath = fullfile(folderName, 'parameters.txt');
discp = fopen(filePath, 'w' );
fprintf(discp, "Experiment Folder: %s\n", rnumb);
fprintf(discp,"Dims: l = %d, p = %d, n = %d, m = %d\n", l, p, n, m);
fprintf(discp,"Trials: %d, Iters: %d\n", num_trials, num_its);
fprintf(discp,"Corr. dist. = %d, %d\n", cor_size(1), cor_size(2));
fprintf(discp,"beta-row = %.5f\n", k_array/m);
fprintf(discp,"beta = %.5f\n",num_corrupt_array/(m*p*n));
fprintf(discp, "q = %.5f\n", q_array);
fprintf(discp, "Wall-clock time (in sec):%.3f\n", tmr);
fclose(discp);

close all

end