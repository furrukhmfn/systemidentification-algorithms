function CompareMeanRMSE()
% CompareMeanRMSE Loads optimization workspaces and plots population mean RMSE
% convergence histories vs iteration across case groups in IEEE publication format.

%% Configuration
dataDir = 'output';
filePattern = 'ALO_Case%d_24.mat';

%% Figure Formatting Defaults
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 9);
set(0, 'DefaultAxesLineWidth', 0.8);
set(0, 'DefaultAxesTickDir', 'in');
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultLineLineWidth', 0.8);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 9);
set(0, 'DefaultLegendFontName', 'Times New Roman');
set(0, 'DefaultLegendFontSize', 7);
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');

%% Fitness History Loader
    function [meanRMSE, maxIt] = loadMeanRMSE(caseNum)
        fname = fullfile(dataDir, sprintf(filePattern, caseNum));
        if ~isfile(fname), error('File not found: %s', fname); end
        S = load(fname, 'globalFitnessMatrix', 'maxIterations', 'iter');
        if ~isfield(S, 'globalFitnessMatrix')
            error('globalFitnessMatrix missing in %s', fname);
        end
        if isfield(S, 'iter'), nIt = S.iter;
        elseif isfield(S, 'maxIterations'), nIt = S.maxIterations;
        else, nIt = size(S.globalFitnessMatrix, 1); end
        maxIt = nIt;
        meanRMSE = zeros(nIt, 1);
        for it = 1:nIt
            meanRMSE(it) = mean(S.globalFitnessMatrix(it, :));
        end
    end

%% Test Groups
groups = {
    'ALO Cases (1,1) - (1,5)',        [11 12 13 14 15];
    'ALO Cases (2,1) - (2,5)',        [21 22 23 24 25];
    'ALO Cases (3,1) - (3,5)',        [31 32 33 34 35];
    'ALO Cases (4,1) - (4,5)',        [41 42 43 44 45];
    'Case (1,1) vs Case (4,1)', [11 41];
    };

colors = lines(10);
styles = {'-', '--', '-.', ':', '--'};

%% Generate Group Convergence Plots
for g = 1:size(groups, 1)
    groupLabel = groups{g, 1};
    caseList   = groups{g, 2};
    nCases     = numel(caseList);

    allMean = cell(nCases, 1);
    allMaxIt = 0;
    for k = 1:nCases
        [m, nIt] = loadMeanRMSE(caseList(k));
        allMean{k} = m;
        allMaxIt = max(allMaxIt, nIt);
    end

    M = NaN(allMaxIt, nCases);
    for k = 1:nCases
        M(1:numel(allMean{k}), k) = allMean{k};
    end

    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.8], 'Color', 'w');
    hold on;
    for k = 1:nCases
        plot(1:allMaxIt, M(:, k), 'LineStyle', styles{mod(k-1,5)+1}, ...
            'Color', colors(k,:), 'LineWidth', 1.2);
    end
    hold off;
    xlabel('Iteration');
    ylabel('Mean RMSE');
    title(groupLabel, 'FontWeight', 'normal', 'FontSize', 9);
    caseLabels = cell(nCases, 1);
    for k = 1:nCases
        caseLabels{k} = sprintf('Case (%d,%d)', floor(caseList(k)/10), mod(caseList(k),10)); 
    end
    legend(caseLabels, ...
        'Location', 'northeast', 'FontSize', 7);
    print(fig, fullfile(dataDir, sprintf('MeanRMSE_%s.png', groupLabel)), '-dpng', '-r600');
    print(fig, fullfile(dataDir, sprintf('MeanRMSE_%s.eps', groupLabel)), '-depsc', '-r600');
end

fprintf('Done. %d figures saved to %s/\n', size(groups,1), dataDir);
end
