function QTRK_plots(tdims, num_corrupt_array, k_array, q_array, num_trials, num_its, cor_size)

%% Code to generate QTRK plots %%

% Create a new folder to the save figure
rnumb = num2str(floor(1000 * rand));
folderName = ['Exp_',num2str(rnumb)] ;
if ~exist(folderName, 'dir') % check it does not exist
    mkdir(folderName);
end

%% Running Algorithm with Different Combination of Parameters

% Number of parameters
param1 = length(num_corrupt_array);
param2 = length(k_array);
param3 = length(q_array);

% Initialize variables to determine y-axis limits
allY = zeros(param1*param2,num_its+1,param3);

% Loop through each combination of parameters to compute aggregated results
for i = 1:param1
    for j = 1:param2
        for l = 1:param3

            % Define current parameters
            num_corrupt = num_corrupt_array(i);
            q = q_array(l);
            k = k_array(j);
    
            % Run QTRK
            errs_matrix = QTRK_trials(tdims, num_corrupt, q, k, cor_size, num_trials, num_its);
    
            % Aggregate results
            % --- >
            % Example: Calculate median errors across trials
            median_errs = median(errs_matrix);
    
            % Collect y-values for later use
            allY((i-1)*param2 + j,:,l) = median_errs;
        end
    end
end

%% Enforcing same y-axis range for all semi-log plots

% Determine global y-axis limits
yMin = min(allY(:));
yMax = max(allY(:));

% Initialize index counter for accessing allY
index = 1;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],[0.2, 0.2, 0.2]};
markers = {'o', '*', '+'};
lineStyles = {'-', '--',':', '-.'};

% Generate subplots of different parameter combinations
for i = 1:param1
    for j = 1:param2

        % Define subplot index
        subplotIndex = (i-1)*param2 + j;

        % Create a subplot
        subplot(param1, param2, subplotIndex);

        individual_fig = figure;
        hold on
        for l = 1:param3
            median_errs = allY(index,:,l);
    
            plot(1:num_its+1, median_errs, 'Color', colors{l},'LineStyle', lineStyles{l}, 'LineWidth', 4);
        end

        % Set the y-axis to a logarithmic scale
        set(gca, 'YScale', 'log');

        % Set font size for tick labels and texts
        set(gca, 'FontSize', 30); %tick labels
        xlabel('Iteration', 'interpreter','latex', FontSize=32);
        ylabel('Relative Error', 'interpreter','latex', FontSize=32);
        legend(['q =', num2str(round(q_array(1),2))],['q =', num2str(round(q_array(2),2))],['q =', num2str(round(q_array(3),2))],['q =', num2str(round(q_array(4),2))],'interpreter','latex', 'Location','southwest',FontSize=32);
        %set(lgnd,'interpreter','latex');
        %set(lgnd,)

        % Set x- and y-axis limits
        ylim([yMin, yMax]);
        xlim([0, num_its]);

        % Save individual subplots
        savefig(individual_fig,fullfile(folderName, sprintf('subplot_%d.fig', subplotIndex)));
        saveas(individual_fig,fullfile(folderName, sprintf('subplot_%d.png', subplotIndex)));
        close(individual_fig);
        hold off
        
        % Update the index for the next entry
        index = index + 1;

    end
end

disp("Experiment Folder:")
disp(rnumb);
close all

end