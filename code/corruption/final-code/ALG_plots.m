function ALG_plots(alg, tdims, num_corrupt_array, k_array, q_array, num_trials, num_its, cor_size, corr_option)

%% Code to generate QTRK or mQTRK plots %%
tic

% Create a new folder to the save figure
rnumb = num2str(floor(1000 * rand));
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
disp("Algorithm: " + alg);
disp("Corr. dist. = " + cor_size);
disp("beta-row = " + k_array/m);
disp("beta = " + num_corrupt_array/(m*p*n));
disp("q = " + q_array);

%% Running Algorithm with Different Combination of Parameters

% Number of parameters
param1 = length(num_corrupt_array);
param2 = length(k_array);
param3 = length(q_array);

% Initialize variables to determine y-axis limits
allY = zeros(param1*param2,num_its+1,param3);

% Loop through each combination of parameters to compute aggregated results
indd = 1;
for i = 1:param1
    % Define current parameters
    num_corrupt = num_corrupt_array(i);

    for j = 1:param2
        % Define current parameters
        k = k_array(j);

        for l = 1:param3
            % Define current parameters
            q = q_array(l);

            disp(i + " out of " + param1 + " and " + j + " out of " + param2 + " and " + l + " out of " + param3)
    
            % Run Algorithm (QTRK or mQTRK)
            errs_matrix = ALG_trials(alg, tdims, num_corrupt, q, k, cor_size, num_trials, num_its);
    
            % Aggregate results
            % --- >
            % Example: Calculate median errors across trials
            median_errs = median(errs_matrix);
    
            % Collect y-values for later use
            allY(indd,:,l) = median_errs;
        end
        indd = indd  + 1;
    end
end

%% Enforcing same y-axis range for all semi-log plots

% Determine global y-axis limits
yMin = min(allY(:));
yMax = max(allY(:));

% Choice of markers, colors, and lines for plotting
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],[0.2, 0.2, 0.2]};
lineStyles = {'-', '--',':', '-.'};
legendLabels = cell(1, param3);

% Initialize index counter for accessing allY
index = 1;

for ll = 1:param3
    legendLabels{ll} = strcat('$q =$ ', num2str(round(q_array(ll),4)));
end

% Generate figures for different parameter combinations
for i = 1:param1
    for j = 1:param2

        individual_fig = figure;
        hold on

        for l = 1:param3
            median_errs = allY(index,:,l);
            plot(1:num_its+1, median_errs, 'Color', colors{l},'LineStyle', lineStyles{l}, 'LineWidth', 6);
        end

        % Set the y-axis to a logarithmic scale
        set(gca, 'YScale', 'log');

        % Set font size for tick labels and texts
        set(gca, 'FontSize', 34); %tick labels
        xlabel('Iteration', 'interpreter','latex', FontSize=38);
        ylabel('Relative Error', 'interpreter','latex', FontSize=38);
        
        if  (i == 1) && (j==1) % show legend only once 
            legend(legendLabels,'Interpreter','latex', 'FontSize', 34);
        end

        % Set x- and y-axis limits
        ylim([yMin, yMax]);
        xlim([0, num_its]);

        % Save individual subfigures
        corr_option = char(corr_option);
        alg = char(alg);
        figFileName = fullfile(folderName, [alg, '_', corr_option, '_exp_', num2str(rnumb), '_subfig_', num2str(num_corrupt_array(i)), '_', num2str(k_array(j)), '.fig']);
        savefig(individual_fig, figFileName);

        % Save individual figures
        set(gcf, 'Position', [100, 100, 600, 400]);  % [left, bottom, width, height]
        pngFileName = fullfile(folderName, [alg, '_', corr_option, '_exp_', num2str(rnumb), '_subfig_', num2str(num_corrupt_array(i)), '_', num2str(k_array(j)), '.png']);
        print(gcf, pngFileName, '-dpng', '-r300');  % Adjust resolution as needed

        close(individual_fig);
        hold off
        
        % Update the index for the next entry
        index = index + 1;

    end
end
tmr = toc;
disp("Wall-clock time (in sec): "  + tmr)

filePath = fullfile(folderName, 'parameters.txt');
discp = fopen(filePath, 'w' );
fprintf(discp, "Experiment Folder: %s\n", rnumb);
fprintf(discp, "Algorithm: %s\n", alg);
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