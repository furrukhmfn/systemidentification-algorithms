function GenerateComparisonGraphics(algorithmNames, outputDir)
% GenerateComparisonGraphics Produces a publication-ready suite of IEEE-format
% comparison figures evaluating all metaheuristic algorithms against the Ho-Kalman baseline.
%
% Figures generated:
%   Fig 1  - Delta-R² Heatmap (Algorithm minus Ho-Kalman)
%   Fig 2  - R² Distribution Box Plots
%   Fig 3  - RMSE Distribution Box Plots
%   Fig 4  - Grouped Mean R² by Population Size
%   Fig 5  - Per-State Mean R² Comparison
%   Fig 6  - Normalised Multi-Metric Radar Profile
%   Fig 7  - Mean Population RMSE Convergence Histories
%   Fig 8  - Parity Plot (Ho-Kalman R² vs Algorithm R²)
%   Fig 9  - Win / Draw / Loss Summary Bar Chart
%   Fig 10 - Mean Algorithm Rank Across All Cases
%   Fig 11 - Absolute R² Performance Heatmap
%   Fig 12 - Per-State Delta-R² Breakdown
%
% Syntax:
%   GenerateComparisonGraphics()
%   GenerateComparisonGraphics({'ALO','GWO','CMAES'})
%   GenerateComparisonGraphics({'ALO','GWO'}, 'output/figs')

%% Setup and Defaults
if nargin < 1 || isempty(algorithmNames)
    algorithmNames = {'ALO','GWO','GHOA','WOA','SSA','SCA','WCA','MFO', ...
                      'CMAES','JADE','SHADE','LSHADE'};
end
if nargin < 2 || isempty(outputDir)
    outputDir = 'output';
end
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

baselineName = 'HOKALMAN';

%% Test Matrix Specification
caseSpecs = [
    11,2000,10000,20; 21,2000,10000,20; 31,2000,10000,20; 41,2000,10000,20;
    12,4000, 8000,20; 22,4000, 8000,20; 32,4000, 8000,20; 42,4000, 8000,20;
    13,6000, 6000,20; 23,6000, 6000,20; 33,6000, 6000,20; 43,6000, 6000,20;
    14,8000, 6000,20; 24,8000, 6000,20; 34,8000, 6000,20; 44,8000, 6000,20;
    15,10000,4000,20; 25,10000,4000,20; 35,10000,4000,20; 45,10000,4000,20;
];
nCases   = size(caseSpecs, 1);
caseNums = caseSpecs(:, 1);

groupLabels = {'N=2000','N=4000','N=6000','N=8000','N=10000'};
groupSize   = 4;

stateNames  = {'u','w','q','\theta','h'};
nStates     = 5;
nAlgos      = numel(algorithmNames);

caseLabels  = arrayfun(@(n) sprintf('(%d,%d)', floor(n/10), mod(n,10)), ...
    caseNums, 'UniformOutput', false);

%% Headless Rendering Configuration
try
    opengl software;
catch
end

%% Figure Formatting Defaults
set(0,'DefaultFigureVisible','off');
set(0,'DefaultFigureColor','w');
set(0,'DefaultAxesColor','w');
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',9);
set(0,'DefaultAxesLineWidth',0.8);
set(0,'DefaultAxesTickDir','in');
set(0,'DefaultAxesBox','on');
set(0,'DefaultLineLineWidth',0.8);
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultTextFontSize',9);
set(0,'DefaultLegendFontName','Times New Roman');
set(0,'DefaultLegendFontSize',8);
set(0,'DefaultTextColor','k');
set(0,'DefaultAxesXColor','k');
set(0,'DefaultAxesYColor','k');

%% Color Palette Setup
nColors = max(nAlgos + 1, 13);
rawCols = [
    0.122 0.471 0.706;  % blue
    0.890 0.102 0.110;  % red
    0.173 0.627 0.173;  % green
    0.584 0.404 0.741;  % purple
    1.000 0.498 0.000;  % orange
    0.651 0.337 0.157;  % brown
    0.969 0.506 0.749;  % pink
    0.500 0.500 0.500;  % grey
    0.737 0.741 0.133;  % olive
    0.090 0.745 0.812;  % cyan
    0.400 0.000 0.600;  % violet
    0.000 0.500 0.000;  % dark green
    0.000 0.000 0.000;  % black (HK)
];
colPad = repmat(rawCols, ceil(nColors/size(rawCols,1)), 1);
algoColors = colPad(1:nAlgos, :);
hkColor    = [0 0 0];

%% Data Collection and Caching
cacheTag  = strjoin(sort(algorithmNames), '_');
cacheFile = fullfile('output', ['_GraphicsCache_' matlab.lang.makeValidName(cacheTag) '.mat']);

metricsAlgo = NaN(nAlgos, nCases, nStates, 6);
metricsHK   = NaN(nCases, nStates, 6);
maxIter     = 20;
convergence = NaN(nAlgos, nCases, maxIter);

if exist(cacheFile, 'file')
    fprintf('\n[CACHE] Loading pre-computed metrics from:\n  %s\n', cacheFile);
    C = load(cacheFile);

    if isequal(C.algorithmNames, algorithmNames) && ...
       isequal(C.caseNums,       caseNums)
        metricsAlgo = C.metricsAlgo;
        metricsHK   = C.metricsHK;
        convergence = C.convergence;
        fprintf('[CACHE] Loaded successfully. Skipping ODE computation.\n\n');
    else
        fprintf('[CACHE] Algorithm list mismatch — recomputing...\n\n');
        cacheFile = '';
    end
end

if ~exist(cacheFile,'file') || ~isequal(exist(cacheFile,'file'),2)
    metricsAlgo = NaN(nAlgos, nCases, nStates, 6);
    metricsHK   = NaN(nCases, nStates, 6);
    convergence = NaN(nAlgos, nCases, maxIter);
end

if all(isnan(metricsHK(:)))
    fprintf('\n========================================================\n');
    fprintf('  GenerateComparisonGraphics — computing metrics...\n');
    fprintf('========================================================\n');

    % Ho-Kalman baseline evaluation
    fprintf('\nPhase 1: Ho-Kalman baseline\n');
    for cIdx = 1:nCases
        f = fullfile('output', sprintf('%s_Case%d_24.mat', baselineName, caseNums(cIdx)));
        if ~exist(f,'file'), fprintf('  [SKIP] HK Case %d not found\n', caseNums(cIdx)); continue; end
        try
            [Xr, Xs] = runODE(f);
            for i = 1:nStates
                metricsHK(cIdx, i, :) = computeMetrics(Xr(i,:)', Xs(:,i));
            end
            fprintf('  HK Case %d OK\n', caseNums(cIdx));
        catch ME
            fprintf('  HK Case %d ERROR: %s\n', caseNums(cIdx), ME.message);
        end
    end

    % Metaheuristic algorithms evaluation
    fprintf('\nPhase 2: Algorithm metrics\n');
    for aIdx = 1:nAlgos
        fprintf('\n  [%d/%d] %s\n', aIdx, nAlgos, algorithmNames{aIdx});
        for cIdx = 1:nCases
            f = fullfile('output', sprintf('%s_Case%d_24.mat', algorithmNames{aIdx}, caseNums(cIdx)));
            if ~exist(f,'file'), continue; end
            try
                [Xr, Xs] = runODE(f);
                for i = 1:nStates
                    metricsAlgo(aIdx, cIdx, i, :) = computeMetrics(Xr(i,:)', Xs(:,i));
                end
                S = load(f, 'globalFitnessMatrix', 'globalrmse', 'iter');
                actualIter = S.iter;
                if isfield(S,'globalrmse') && numel(S.globalrmse) >= actualIter
                    cv = S.globalrmse(1:actualIter);
                else
                    cv = mean(S.globalFitnessMatrix(1:actualIter,:), 2);
                end
                convergence(aIdx, cIdx, 1:actualIter) = cv(:);
            catch ME
                fprintf('    Case %d ERROR: %s\n', caseNums(cIdx), ME.message);
            end
        end
    end

    cacheFile = fullfile('output', ['_GraphicsCache_' matlab.lang.makeValidName(cacheTag) '.mat']);
    fprintf('\n[CACHE] Saving metrics to cache:\n  %s\n', cacheFile);
    save(cacheFile, 'metricsAlgo', 'metricsHK', 'convergence', 'algorithmNames', 'caseNums', '-v7.3');
    fprintf('[CACHE] Saved.\n\n');
end

%% Derived Metrics Calculation
r2Algo  = squeeze(mean(metricsAlgo(:,:,:,5), 3, 'omitnan'));
fitAlgo = squeeze(mean(metricsAlgo(:,:,:,6), 3, 'omitnan'));
rmseAlgo= squeeze(mean(metricsAlgo(:,:,:,1), 3, 'omitnan'));

r2HK    = squeeze(mean(metricsHK(:,:,5), 2, 'omitnan'));
fitHK   = squeeze(mean(metricsHK(:,:,6), 2, 'omitnan'));
rmseHK  = squeeze(mean(metricsHK(:,:,1), 2, 'omitnan'));

deltaR2  = r2Algo   - r2HK';
deltaFIT = fitAlgo  - fitHK';
deltaRMSE= rmseAlgo - rmseHK';

fprintf('\nAll metrics ready. Generating figures...\n\n');

%% Fig 1: Delta-R2 Heatmap
fig = figure('Units','inches','Position',[0 0 7.16 4.0],'Color','w');
ax  = axes(fig);

clim = max(abs(deltaR2(:)), [], 'omitnan');
clim = max(clim, 0.01);

imagesc(ax, deltaR2);
set(ax, 'XTick', 1:nCases, 'XTickLabel', caseLabels, ...
        'XTickLabelRotation', 45, ...
        'YTick', 1:nAlgos, 'YTickLabel', algorithmNames, ...
        'FontName','Times New Roman','FontSize',8, ...
        'CLim', [-clim clim]);

n  = 128;
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); ...
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];
colormap(ax, cmap);
cb = colorbar(ax);
cb.Label.String = '$\Delta R^2$ (Algorithm $-$ Ho-Kalman)';
cb.Label.Interpreter = 'latex';
cb.Label.FontName = 'Times New Roman';  cb.Label.FontSize = 9;

for aIdx = 1:nAlgos
    for cIdx = 1:nCases
        v = deltaR2(aIdx, cIdx);
        if ~isnan(v)
            tc = 'k';
            if abs(v) > 0.5*clim, tc = 'w'; end
            text(ax, cIdx, aIdx, sprintf('%.2f', v), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',6,'Color',tc,'FontName','Times New Roman');
        end
    end
end

hold(ax,'on');
for g = 1:4
    xline(ax, g*4 + 0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end
xlabel(ax,'Case');
ylabel(ax,'Algorithm');
title(ax,'$\Delta R^2$ vs Ho-Kalman Baseline (mean over 5 states)', ...
    'Interpreter','latex','FontWeight','normal','FontSize',9);
saveFig(fig, outputDir, 'Comparison_Heatmap_DeltaR2');

%% Fig 2: R2 Distribution Box Plot
fig = figure('Units','inches','Position',[0 0 7.16 3.5],'Color','w');
ax  = axes(fig);
hold(ax,'on');

for aIdx = 1:nAlgos
    drawBoxPlot(ax, aIdx, r2Algo(aIdx,:)', algoColors(aIdx,:), 0.35);
end

hkMeanAll = mean(r2HK, 'omitnan');
yline(ax, hkMeanAll, '--k', 'LineWidth', 1.2);
text(ax, 0.6, hkMeanAll, sprintf(' HK mean = %.3f', hkMeanAll), ...
    'VerticalAlignment','bottom','FontName','Times New Roman','FontSize',7);

set(ax,'XTick',1:nAlgos,'XTickLabel',algorithmNames,'XTickLabelRotation',30);
xlim(ax,[0.5 nAlgos+0.5]);
ylab = max(0, min(r2Algo(:),[],'omitnan') - 0.05);
ylim(ax,[ylab, 1.0]);
xlabel(ax,'Algorithm');
ylabel(ax,'R^2 (mean over 5 states)');
title(ax,'R^2 Distribution Across 20 Cases', ...
    'FontWeight','normal','FontSize',9);
grid(ax,'on'); ax.GridAlpha = 0.2;
hold(ax,'off');
saveFig(fig, outputDir, 'Comparison_BoxPlot_R2');

%% Fig 3: RMSE Distribution Box Plot
fig = figure('Units','inches','Position',[0 0 7.16 3.5],'Color','w');
ax  = axes(fig);
hold(ax,'on');

for aIdx = 1:nAlgos
    drawBoxPlot(ax, aIdx, rmseAlgo(aIdx,:)', algoColors(aIdx,:), 0.35);
end

hkRmseMean = mean(rmseHK, 'omitnan');
yline(ax, hkRmseMean, '--k', 'LineWidth', 1.2);
text(ax, 0.6, hkRmseMean, sprintf(' HK mean = %.3f', hkRmseMean), ...
    'VerticalAlignment','bottom','FontName','Times New Roman','FontSize',7);

set(ax,'XTick',1:nAlgos,'XTickLabel',algorithmNames,'XTickLabelRotation',30);
xlim(ax,[0.5 nAlgos+0.5]);
xlabel(ax,'Algorithm'); ylabel(ax,'RMSE (mean over 5 states)');
title(ax,'RMSE Distribution Across 20 Cases', ...
    'FontWeight','normal','FontSize',9);
grid(ax,'on'); ax.GridAlpha = 0.2;
hold(ax,'off');
saveFig(fig, outputDir, 'Comparison_BoxPlot_RMSE');

%% Fig 4: Grouped Bar Chart by Population Size
fig = figure('Units','inches','Position',[0 0 7.16 4.5],'Color','w');

nGroups = 5;
for gIdx = 1:nGroups
    ax = subplot(1, nGroups, gIdx);
    cIdx1 = (gIdx-1)*4 + 1;
    cIdx2 = gIdx*4;
    cIdxRange = cIdx1:cIdx2;

    r2group = mean(r2Algo(:, cIdxRange), 2, 'omitnan');
    r2hkGrp = mean(r2HK(cIdxRange),     'omitnan');

    [r2sort, sortI] = sort(r2group, 'descend');
    namSort = algorithmNames(sortI);
    colSort = algoColors(sortI, :);

    bh = barh(ax, r2sort, 0.7, 'FaceColor','flat');
    bh.CData = colSort;
    xline(ax, r2hkGrp, '--k', 'LineWidth', 1.0);

    set(ax, 'YTick', 1:nAlgos, 'YTickLabel', namSort, ...
        'FontName','Times New Roman','FontSize',7);
    xlabel(ax, 'Mean R^2','FontSize',8);
    title(ax, groupLabels{gIdx}, 'FontWeight','normal','FontSize',8);
    xlim(ax, [max(0, min(r2sort)-0.1), 1.0]);
    grid(ax,'on'); ax.GridAlpha = 0.2;
end

sgtitle('Mean R^2 per Population Group (Ho-Kalman dashed)', ...
    'FontWeight','normal','FontSize',9);
saveFig(fig, outputDir, 'Comparison_GroupedBar_R2byGroup');

%% Fig 5: Per-State Mean R2 Bar Chart
fig = figure('Units','inches','Position',[0 0 7.16 4.0],'Color','w');
ax  = axes(fig);

r2byState = zeros(nStates, nAlgos);
for i = 1:nStates
    for aIdx = 1:nAlgos
        r2byState(i, aIdx) = mean(squeeze(metricsAlgo(aIdx,:,i,5)), 'omitnan');
    end
end

r2hkByState = zeros(nStates, 1);
for i = 1:nStates
    r2hkByState(i) = mean(metricsHK(:,i,5), 'omitnan');
end

bh = bar(ax, r2byState);
for aIdx = 1:nAlgos
    bh(aIdx).FaceColor = algoColors(aIdx,:);
    bh(aIdx).EdgeColor = 'none';
end

hold(ax,'on');
for i = 1:nStates
    plot(ax, i, r2hkByState(i), 'dk', 'MarkerSize', 8, ...
        'MarkerFaceColor','k','LineWidth',1.0);
end

set(ax, 'XTick', 1:nStates, ...
    'XTickLabel', {'$\dot{u}$','$\dot{w}$','$\dot{q}$','$\dot{\theta}$','$\dot{h}$'}, ...
    'TickLabelInterpreter','latex');
ylabel(ax,'Mean R^2 (over 20 cases)');
xlabel(ax,'State');
title(ax,'Per-State Mean R^2 by Algorithm (HK shown as filled diamond)', ...
    'FontWeight','normal','FontSize',9);
ylim(ax,[0 1]);
grid(ax,'on'); ax.GridAlpha = 0.2;

lgdStr  = [algorithmNames, {'Ho-Kalman'}];
hkLgd   = plot(ax, nan, nan, 'dk', 'MarkerFaceColor','k','MarkerSize',8);
lgdObjs = [bh(:); hkLgd];
legend(ax, lgdObjs, lgdStr, 'Location','southeast', 'FontSize',7, ...
    'NumColumns', ceil((nAlgos+1)/3));
hold(ax,'off');
saveFig(fig, outputDir, 'Comparison_PerState_R2');

%% Fig 6: Multi-Metric Radar Profile
fig = figure('Units','inches','Position',[0 0 5.0 5.0],'Color','w');
ax  = axes(fig,'Visible','off');

radarMetrics = {'Mean R^2', 'Mean FIT%', '1/RMSE', '1/MAE', '1/NRMSE'};
nRad         = numel(radarMetrics);

rawVals = zeros(nAlgos+1, nRad);
for aIdx = 1:nAlgos
    rawVals(aIdx, 1) = mean(r2Algo(aIdx,:),   'omitnan');
    rawVals(aIdx, 2) = mean(fitAlgo(aIdx,:),   'omitnan') / 100;
    rawVals(aIdx, 3) = 1 / max(mean(rmseAlgo(aIdx,:), 'omitnan'), eps);
    rawVals(aIdx, 4) = 1 / max(mean(squeeze(mean(metricsAlgo(aIdx,:,:,2),3,'omitnan')), 'omitnan'), eps);
    rawVals(aIdx, 5) = 1 / max(mean(squeeze(mean(metricsAlgo(aIdx,:,:,3),3,'omitnan')), 'omitnan'), eps);
end

rawVals(end, 1) = mean(r2HK,   'omitnan');
rawVals(end, 2) = mean(fitHK,  'omitnan') / 100;
rawVals(end, 3) = 1 / max(mean(rmseHK, 'omitnan'), eps);
rawVals(end, 4) = 1 / max(mean(squeeze(mean(metricsHK(:,:,2),2,'omitnan')), 'omitnan'), eps);
rawVals(end, 5) = 1 / max(mean(squeeze(mean(metricsHK(:,:,3),2,'omitnan')), 'omitnan'), eps);

vMin  = min(rawVals, [], 1);
vMax  = max(rawVals, [], 1);
vRng  = max(vMax - vMin, eps);
normVals = (rawVals - vMin) ./ vRng;

angles  = linspace(0, 2*pi, nRad+1);
angles  = angles(1:end-1);
axX     = cos(angles - pi/2);
axY     = sin(angles - pi/2);

hold(ax,'on'); ax.Visible = 'off';
for rr = 0.25:0.25:1.0
    th = linspace(0,2*pi,200);
    plot(ax, rr*cos(th), rr*sin(th), ':', 'Color',[0.8 0.8 0.8], 'LineWidth',0.5);
end

for k = 1:nRad
    plot(ax, [0 axX(k)], [0 axY(k)], 'k-', 'LineWidth', 0.5);
    text(ax, 1.15*axX(k), 1.15*axY(k), radarMetrics{k}, ...
        'HorizontalAlignment','center','FontName','Times New Roman','FontSize',8);
end

lgdHandles = zeros(nAlgos+1, 1);
for aIdx = 1:nAlgos+1
    v   = normVals(aIdx, :);
    px  = [v .* axX, v(1)*axX(1)];
    py  = [v .* axY, v(1)*axY(1)];
    if aIdx <= nAlgos
        col = algoColors(aIdx,:);
        lw  = 0.8;
        ls  = '-';
    else
        col = hkColor;
        lw  = 1.5;
        ls  = '--';
    end
    h = plot(ax, px, py, ls, 'Color', col, 'LineWidth', lw);
    lgdHandles(aIdx) = h;
end
axis(ax, 'equal'); axis(ax, 'off');
xlim(ax,[-1.4 1.4]); ylim(ax,[-1.4 1.4]);

lgdStr = [algorithmNames, {'Ho-Kalman'}];
legend(ax, lgdHandles, lgdStr, 'Location','southoutside', ...
    'NumColumns', ceil((nAlgos+1)/2), 'FontSize',7);
title(ax,'Normalised Multi-Metric Performance Profile', ...
    'FontWeight','normal','FontSize',9,'Visible','on');
saveFig(fig, outputDir, 'Comparison_Radar_MultiMetric');

%% Fig 7: Population Mean Convergence Curves
fig = figure('Units','inches','Position',[0 0 3.5 3.0],'Color','w');
ax  = axes(fig);

hold(ax,'on');
for aIdx = 1:nAlgos
    cv  = squeeze(convergence(aIdx,:,:));
    cvM = mean(cv, 1, 'omitnan');
    validIter = find(~isnan(cvM));
    if isempty(validIter), continue; end
    plot(ax, validIter, cvM(validIter), '-', ...
        'Color', algoColors(aIdx,:), 'LineWidth', 0.9);
end
hold(ax,'off');

legend(ax, algorithmNames, 'Location','northeast', 'FontSize',7, ...
    'NumColumns', ceil(nAlgos/4));
xlabel(ax,'Iteration');
ylabel(ax,'Mean Population RMSE');
title(ax,'Convergence Curves (mean over 20 cases)', ...
    'FontWeight','normal','FontSize',9);
grid(ax,'on'); ax.GridAlpha = 0.2;
set(ax,'YScale','log');
saveFig(fig, outputDir, 'Comparison_ConvergenceCurves');

%% Fig 8: Parity Scatter Plot
fig = figure('Units','inches','Position',[0 0 3.5 3.5],'Color','w');
ax  = axes(fig);

markerTypes = {'o','s','^','v','d','p','h','+','x','*','<','>','.'};
hold(ax,'on');
lgdHandles2 = zeros(nAlgos,1);
for aIdx = 1:nAlgos
    mk = markerTypes{mod(aIdx-1, numel(markerTypes))+1};
    h = scatter(ax, r2HK, r2Algo(aIdx,:)', 18, ...
        'Marker', mk, 'MarkerEdgeColor', algoColors(aIdx,:), ...
        'MarkerFaceColor', algoColors(aIdx,:), ...
        'MarkerFaceAlpha', 0.6);
    lgdHandles2(aIdx) = h;
end

allR2 = [r2HK; r2Algo(:)];
mn = min(allR2(:), [], 'omitnan') - 0.02;
mx = max(allR2(:), [], 'omitnan') + 0.02;
plot(ax, [mn mx], [mn mx], 'k--', 'LineWidth', 1.0);
text(ax, (mn+mx)/2 + 0.02, (mn+mx)/2 + 0.02, 'HK = Algo', ...
    'FontSize',7,'FontName','Times New Roman','Rotation',45);
hold(ax,'off');

xlabel(ax,'Ho-Kalman R^2');
ylabel(ax,'Algorithm R^2');
title(ax,'Parity Plot: Algorithm vs Ho-Kalman R^2', ...
    'FontWeight','normal','FontSize',9);
legend(ax, lgdHandles2, algorithmNames, 'Location','southeast', ...
    'FontSize',7, 'NumColumns', ceil(nAlgos/4));
axis(ax,'square');
xlim(ax,[mn mx]); ylim(ax,[mn mx]);
grid(ax,'on'); ax.GridAlpha = 0.2;
saveFig(fig, outputDir, 'Comparison_Parity_R2');

%% Fig 9: Win/Draw/Loss Stacked Bar Chart
fig = figure('Units','inches','Position',[0 0 3.5 4.0],'Color','w');
ax  = axes(fig);

wdl = zeros(nAlgos, 3);
for aIdx = 1:nAlgos
    for cIdx = 1:nCases
        d = r2Algo(aIdx,cIdx) - r2HK(cIdx);
        if isnan(d), continue; end
        if d >  0.001, wdl(aIdx,1) = wdl(aIdx,1)+1;
        elseif d < -0.001, wdl(aIdx,3) = wdl(aIdx,3)+1;
        else,          wdl(aIdx,2) = wdl(aIdx,2)+1;
        end
    end
end

[~,si] = sort(wdl(:,1), 'descend');
wdlS   = wdl(si,:);
namS   = algorithmNames(si);

bh = barh(ax, wdlS, 'stacked');
bh(1).FaceColor = [0.173 0.627 0.173];
bh(2).FaceColor = [0.750 0.750 0.750];
bh(3).FaceColor = [0.890 0.102 0.110];

set(ax,'YTick',1:nAlgos,'YTickLabel',namS, ...
    'FontName','Times New Roman','FontSize',8);
xlabel(ax,'Number of Cases (out of 20)');
title(ax,'Win / Draw / Loss vs Ho-Kalman', ...
    'FontWeight','normal','FontSize',9);
legend(ax,bh,{'Win','Draw','Loss'},'Location','southeast','FontSize',8);
xline(ax,10,'--k','LineWidth',0.8);
xlim(ax,[0 nCases]);
grid(ax,'on'); ax.GridAlpha = 0.2; ax.XGrid='on'; ax.YGrid='off';
saveFig(fig, outputDir, 'Comparison_WinDrawLoss');

%% Fig 10: Mean Rank Across Test Cases
fig = figure('Units','inches','Position',[0 0 3.5 3.5],'Color','w');
ax  = axes(fig);

r2All  = [r2Algo; r2HK'];
nameAll= [algorithmNames, {baselineName}];
nAll   = nAlgos + 1;

ranks  = zeros(nAll, nCases);
for cIdx = 1:nCases
    [~,sortI] = sort(r2All(:,cIdx), 'descend');
    r = zeros(nAll,1);
    r(sortI) = 1:nAll;
    ranks(:,cIdx) = r;
end
meanRank = mean(ranks, 2, 'omitnan');
stdRank  = std(ranks,  0, 2, 'omitnan');

[mrS, si] = sort(meanRank, 'ascend');
stdS      = stdRank(si);
namRS     = nameAll(si);

colRS = [algoColors; hkColor];
colRS = colRS(si,:);

barh(ax, mrS, 0.65, 'FaceColor','flat', 'CData', colRS);
hold(ax,'on');
errorbar(ax, mrS, 1:nAll, stdS, 'horizontal', 'k.', 'LineWidth',0.8, 'CapSize',4);
xline(ax, mean(1:nAll), '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
hold(ax,'off');

set(ax,'YTick',1:nAll,'YTickLabel',namRS, ...
    'FontName','Times New Roman','FontSize',8);
xlabel(ax,'Mean Rank (1 = best)');
title(ax,'Algorithm Ranking by R^2 (lower is better)', ...
    'FontWeight','normal','FontSize',9);
grid(ax,'on'); ax.GridAlpha=0.2; ax.XGrid='on'; ax.YGrid='off';
xlim(ax,[0 nAll+0.5]);
saveFig(fig, outputDir, 'Comparison_MeanRank');

%% Fig 11: Absolute R2 Heatmap
fig = figure('Units','inches','Position',[0 0 7.16 4.5],'Color','w');
ax  = axes(fig);

r2plotData = [r2All(1:nAlgos,:); r2HK'];
r2plotNames= [algorithmNames, {[baselineName ' (Baseline)']}];
nRowsFig11 = nAlgos + 1;

imagesc(ax, r2plotData);
colormap(ax, parula(256));
cb = colorbar(ax);
cb.Label.String = 'R^2';
cb.Label.FontName = 'Times New Roman'; cb.Label.FontSize = 9;
caxis(ax, [max(0, min(r2plotData(:), [], 'omitnan') - 0.05), 1.0]);

set(ax,'XTick',1:nCases,'XTickLabel',caseLabels, ...
    'XTickLabelRotation',45, ...
    'YTick',1:nRowsFig11,'YTickLabel',r2plotNames, ...
    'FontName','Times New Roman','FontSize',8);

for r = 1:nRowsFig11
    for c = 1:nCases
        v = r2plotData(r,c);
        if ~isnan(v)
            tc = 'k';
            if v < 0.3, tc='w'; end
            text(ax,c,r,sprintf('%.2f',v), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',6,'Color',tc,'FontName','Times New Roman');
        end
    end
end

hold(ax,'on');
plot(ax,[0.5 nCases+0.5], [nAlgos+0.5 nAlgos+0.5], 'w-','LineWidth',1.5);
for g=1:4
    xline(ax, g*4+0.5,'--','Color',[0.8 0.8 0.8],'LineWidth',0.5);
end
hold(ax,'off');

xlabel(ax,'Case'); ylabel(ax,'Algorithm');
title(ax,'R^2 Heatmap — All Algorithms Including Ho-Kalman Baseline', ...
    'FontWeight','normal','FontSize',9);
saveFig(fig, outputDir, 'Comparison_Heatmap_AbsoluteR2');

%% Fig 12: Per-State Delta-R2 Breakdown
fig = figure('Units','inches','Position',[0 0 7.16 8.0],'Color','w');
stateFullNames = {'$u$ (m/s)','$w$ (m/s)','$q$ (rad/s)','$\theta$ (rad)','$h$ (m)'};

for sIdx = 1:nStates
    ax = subplot(3, 2, sIdx);

    dR2state = zeros(nAlgos, 1);
    for aIdx = 1:nAlgos
        dR2state(aIdx) = mean(squeeze(metricsAlgo(aIdx,:,sIdx,5)) - metricsHK(:,sIdx,5)', 'omitnan');
    end

    [dSort, si] = sort(dR2state, 'descend');
    namSort      = algorithmNames(si);
    colSort      = algoColors(si,:);

    bh = barh(ax, dSort, 0.7, 'FaceColor','flat');
    bh.CData = colSort;
    xline(ax, 0, 'k-', 'LineWidth', 1.0);

    posIdx = dSort > 0.001;
    negIdx = dSort < -0.001;
    if any(posIdx)
        bh.CData(posIdx,:)  = repmat([0.173 0.627 0.173], sum(posIdx), 1);
    end
    if any(negIdx)
        bh.CData(negIdx,:)  = repmat([0.890 0.102 0.110], sum(negIdx), 1);
    end

    set(ax,'YTick',1:nAlgos,'YTickLabel',namSort, ...
        'FontName','Times New Roman','FontSize',7);
    xlabel(ax,'$\Delta R^2$ (Algo $-$ HK)','Interpreter','latex','FontSize',8);
    title(ax, stateFullNames{sIdx}, 'Interpreter','latex', ...
        'FontWeight','normal','FontSize',9);
    grid(ax,'on'); ax.GridAlpha=0.2; ax.XGrid='on'; ax.YGrid='off';
    ax.XAxisLocation='bottom';
end

ax = subplot(3,2,6); axis(ax,'off');
p1 = patch(ax,[0 1 1 0],[0 0 1 1],[0.173 0.627 0.173],'EdgeColor','none');
p2 = patch(ax,[0 1 1 0],[0 0 1 1],[0.890 0.102 0.110],'EdgeColor','none');
legend(ax,[p1,p2],{'Better than HK','Worse than HK'}, ...
    'Location','best','FontSize',9);

sgtitle('$\Delta R^2$ vs Ho-Kalman per State (green = algo wins, red = HK wins)', ...
    'Interpreter','latex','FontWeight','normal','FontSize',9);
saveFig(fig, outputDir, 'Comparison_PerState_DeltaR2');

fprintf('\n========================================================\n');
fprintf('  All figures saved to: %s/\n', outputDir);
fprintf('========================================================\n');
end

%% Local Helper Functions

function saveFig(fig, outDir, baseName)
% Saves figure as 600-dpi raster PNG and vector/image PDF
    pngFile = fullfile(outDir, [baseName '.png']);
    pdfFile = fullfile(outDir, [baseName '.pdf']);
    savedFmts = {};

    pngOK = false;
    try
        exportgraphics(fig, pngFile, 'Resolution', 600);
        pngOK = true;
    catch
    end
    if ~pngOK
        try
            set(fig, 'PaperPositionMode', 'auto');
            print(fig, '-dpng', '-r600', pngFile);
            pngOK = true;
        catch ME2
            fprintf('    [WARN] PNG save failed: %s\n', ME2.message);
        end
    end
    if pngOK, savedFmts{end+1} = '.png'; end

    pdfOK = false;
    try
        exportgraphics(fig, pdfFile, 'ContentType', 'image', 'Resolution', 600);
        pdfOK = true;
    catch
    end
    if ~pdfOK
        try
            set(fig, 'PaperPositionMode', 'auto');
            print(fig, '-dpdf', pdfFile);
            pdfOK = true;
        catch ME2
            fprintf('    [WARN] PDF save failed: %s\n', ME2.message);
        end
    end
    if pdfOK, savedFmts{end+1} = '.pdf'; end

    fprintf('  Saved: %s  (%s)\n', baseName, strjoin(savedFmts, ' + '));
    close(fig);
end

function [X_real, X_sim] = runODE(matFile)
% Solves longitudinal state equations over flight data window using ode15s
    S             = load(matFile);
    Vb            = S.Vb;
    pqrD          = S.pqr;
    phi_theta_psi = S.phi_theta_psi;
    Xe            = S.Xe;
    elevator      = S.elevator;
    flapPos       = S.flapPos;
    flapNeg       = S.flapNeg;
    flapDiff      = S.flapDiff;
    globalBest    = S.globalBest;
    staticData    = S.staticData;
    try, margin = S.inputData.ValidationMargin; catch, margin = 0; end

    N = size(Vb,1);  dS = 1+margin;  dE = N-margin;
    if isfield(staticData,'time_all'), time = staticData.time_all;
    else, time = (0:N-1)'; time = time(dS:dE); end

    u_d  = Vb(dS:dE,1);   w_d  = Vb(dS:dE,3);
    q_d  = pqrD(dS:dE,2); th_d = phi_theta_psi(2,dS:dE)';
    h_d  = -Xe(dS:dE,3);
    de   = elevator(dS:dE,1);
    fp   = flapPos(dS:dE,1);
    fn   = flapNeg(dS:dE,1);
    fd   = flapDiff(dS:dE,1);

    X_real = [u_d,w_d,q_d,th_d,h_d]';
    [A,B]  = formatParameters24(globalBest, mean(u_d));
    odefun = @(t,x) ode_tv(t,x,A,B,time,u_d,de,fp,fn,fd,0);
    opts   = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [ts,Xr]= ode15s(odefun, time, X_real(:,1), opts);
    X_sim  = interp1(ts, Xr, time, 'linear','extrap');
end

function dx = ode_tv(t,x,A,B,time,u,elev,fp,fn,fd,thrust)
% Evaluates time-varying longitudinal system dynamics
    de=interp1(time,elev,t,'linear','extrap');
    fpp=interp1(time,fp,t,'linear','extrap');
    fn_=interp1(time,fn,t,'linear','extrap');
    fd_=interp1(time,fd,t,'linear','extrap');
    ut=interp1(time,u,t,'linear','extrap');
    Atv=A; Atv(5,4)=ut;
    dx=Atv*x+B*[de;thrust;fpp;fn_;fd_];
end

function m = computeMetrics(y, yhat)
% Returns [RMSE, MAE, NRMSE, MAPE, R2, FIT]
    y=y(:); yhat=yhat(:);
    v=isfinite(y)&isfinite(yhat); y=y(v); yhat=yhat(v);
    if isempty(y), m=NaN(1,6); return; end
    err=y-yhat;
    RMSE=sqrt(mean(err.^2)); MAE=mean(abs(err));
    yr=max(y)-min(y); NRMSE=NaN; if yr~=0, NRMSE=RMSE/yr; end
    MAPE=100*mean(abs(err)./max(abs(y),eps*max(abs(y))));
    SS_res=sum(err.^2); SS_tot=sum((y-mean(y)).^2);
    R2=NaN; FIT=NaN;
    if SS_tot~=0, R2=1-SS_res/SS_tot; FIT=100*(1-sqrt(SS_res/SS_tot)); end
    m=[RMSE,MAE,NRMSE,MAPE,R2,FIT];
end

function drawBoxPlot(ax, xPos, data, col, halfWidth)
% Custom base-MATLAB box-and-whisker plot rendering
    data = data(isfinite(data));
    if isempty(data), return; end

    q1  = prctile(data, 25);
    q2  = median(data);
    q3  = prctile(data, 75);
    iqr_ = q3 - q1;
    wLow = max(data(data >= q1 - 1.5*iqr_));
    wHig = min(data(data <= q3 + 1.5*iqr_));
    if isempty(wLow), wLow = q1; end
    if isempty(wHig), wHig = q3; end
    outliers = data(data < wLow | data > wHig);

    x1 = xPos - halfWidth;
    x2 = xPos + halfWidth;

    patch(ax, [x1 x2 x2 x1 x1], [q1 q1 q3 q3 q1], col, ...
        'FaceAlpha', 0.35, 'EdgeColor', col, 'LineWidth', 0.8);

    plot(ax, [x1 x2], [q2 q2], '-', 'Color', col, 'LineWidth', 1.8);

    mn = mean(data);
    plot(ax, xPos, mn, 'd', 'MarkerSize', 5, ...
        'MarkerEdgeColor', col, 'MarkerFaceColor', 'w', 'LineWidth', 0.8);

    plot(ax, [xPos xPos], [q1 wLow], '-', 'Color', col, 'LineWidth', 0.8);
    plot(ax, [xPos xPos], [q3 wHig], '-', 'Color', col, 'LineWidth', 0.8);
    plot(ax, [x1+0.05 x2-0.05], [wLow wLow], '-', 'Color', col, 'LineWidth', 0.8);
    plot(ax, [x1+0.05 x2-0.05], [wHig wHig], '-', 'Color', col, 'LineWidth', 0.8);

    if ~isempty(outliers)
        plot(ax, xPos*ones(size(outliers)), outliers, '+', ...
            'Color', col, 'MarkerSize', 5, 'LineWidth', 0.8);
    end
end
