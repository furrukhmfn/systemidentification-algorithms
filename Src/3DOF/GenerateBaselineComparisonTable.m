function GenerateBaselineComparisonTable(algorithmNames, outputExcel)
% GenerateBaselineComparisonTable Evaluates and compares all metaheuristic
% optimization algorithms against the Ho-Kalman (HOKALMAN) algebraic baseline
% across the complete 20-case test matrix for longitudinal aircraft dynamics.
%
% Methodology:
%   1. Loads identified model matrices for Ho-Kalman and metaheuristic algorithms.
%   2. Runs ODE15s state simulations over matching flight test validation windows.
%   3. Computes per-state metrics: RMSE, MAE, NRMSE, MAPE, R², and FIT%.
%   4. Computes delta (Algorithm - Baseline) and ratio (Algorithm / Baseline) metrics.
%   5. Compiles cross-algorithm comparison matrices, per-case rankings, and win/loss records.
%
% Outputs:
%   Excel report written to output/<outputExcel> containing individual algorithm sheets,
%   cross-comparison performance matrices, rankings, and per-state breakdowns.
%
% Syntax:
%   GenerateBaselineComparisonTable()
%   GenerateBaselineComparisonTable({'ALO','GWO','CMAES'})
%   GenerateBaselineComparisonTable({'ALO','GWO'}, 'MyComparison.xlsx')

%% Setup and Defaults
if nargin < 1 || isempty(algorithmNames)
    algorithmNames = {'ALO','GWO','GHOA','WOA','SSA','SCA','WCA','MFO', ...
                      'CMAES','JADE','SHADE','LSHADE'};
end
if nargin < 2 || isempty(outputExcel)
    outputExcel = 'BaselineComparisonTable_vs_HOKALMAN.xlsx';
end

baselineName = 'HOKALMAN';

%% Test Matrix Specification
caseSpecs = [
    11, 2000, 10000, 20;
    21, 2000, 10000, 20;
    31, 2000, 10000, 20;
    41, 2000, 10000, 20;
    12, 4000,  8000, 20;
    22, 4000,  8000, 20;
    32, 4000,  8000, 20;
    42, 4000,  8000, 20;
    13, 6000,  6000, 20;
    23, 6000,  6000, 20;
    33, 6000,  6000, 20;
    43, 6000,  6000, 20;
    14, 8000,  6000, 20;
    24, 8000,  6000, 20;
    34, 8000,  6000, 20;
    44, 8000,  6000, 20;
    15,10000,  4000, 20;
    25,10000,  4000, 20;
    35,10000,  4000, 20;
    45,10000,  4000, 20;
];
nCases   = size(caseSpecs, 1);
caseNums = caseSpecs(:, 1);

%% State and Metric Definitions
stateNames  = {'u (m/s)', 'w (m/s)', 'q (rad/s)', 'theta (rad)', 'h (m)'};
nStates     = 5;
metricNames = {'RMSE', 'MAE', 'NRMSE', 'MAPE_%', 'R2', 'FIT_%'};
nMetrics    = numel(metricNames);
nAlgos      = numel(algorithmNames);

higherBetter = [false false false false true true];

%% Pre-allocation
metricsAlgo = NaN(nAlgos, nCases, nStates, nMetrics);
metricsHK   = NaN(nCases, nStates, nMetrics);

%% Execution Banner
fprintf('\n========================================================\n');
fprintf('  GenerateBaselineComparisonTable\n');
fprintf('  Baseline : %s\n', baselineName);
fprintf('  Algorithms compared : %d\n', nAlgos);
fprintf('  Cases per algorithm : %d\n', nCases);
fprintf('========================================================\n');

%% Phase 1: Ho-Kalman Baseline Evaluation
fprintf('\n--- Phase 1: Computing Ho-Kalman baseline metrics ---\n');

for cIdx = 1:nCases
    caseNum = caseNums(cIdx);
    hkFile  = fullfile('output', sprintf('%s_Case%d_24.mat', baselineName, caseNum));

    if ~exist(hkFile, 'file')
        fprintf('  [SKIP HK] %s not found.\n', hkFile);
        continue
    end

    try
        [X_real, X_sim, time] = runODE(hkFile);
        for i = 1:nStates
            metricsHK(cIdx, i, :) = computeMetrics(X_real(i,:)', X_sim(:,i));
        end
        fprintf('  HK Case %d - OK\n', caseNum);
    catch ME
        fprintf('  HK Case %d - ERROR: %s\n', caseNum, ME.message);
    end
end

%% Phase 2: Metaheuristic Algorithm Evaluation
fprintf('\n--- Phase 2: Computing algorithm metrics ---\n');

for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    fprintf('\n[%d/%d] %s\n', aIdx, nAlgos, algoName);

    for cIdx = 1:nCases
        caseNum  = caseNums(cIdx);
        algoFile = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, caseNum));

        if ~exist(algoFile, 'file')
            fprintf('  [SKIP] %s not found.\n', algoFile);
            continue
        end

        try
            [X_real, X_sim, ~] = runODE(algoFile);
            for i = 1:nStates
                metricsAlgo(aIdx, cIdx, i, :) = computeMetrics(X_real(i,:)', X_sim(:,i));
            end
            fprintf('  Case %d (%d/%d) - OK\n', caseNum, cIdx, nCases);
        catch ME
            fprintf('  Case %d (%d/%d) - ERROR: %s\n', caseNum, cIdx, nCases, ME.message);
        end
    end
end

%% Metric Delta and Ratio Calculations
hkExpanded = repmat(reshape(metricsHK, [1, nCases, nStates, nMetrics]), [nAlgos 1 1 1]);
deltaMetrics = metricsAlgo - hkExpanded;
ratioMetrics = metricsAlgo ./ (hkExpanded + eps);

%% Excel Export Setup
if ~exist('output', 'dir'), mkdir('output'); end
excelOut = fullfile('output', outputExcel);
if exist(excelOut, 'file'), delete(excelOut); end
fprintf('\nWriting results to: %s\n', excelOut);

caseLabels = arrayfun(@(n) sprintf('(%d,%d)', floor(n/10), mod(n,10)), ...
    caseNums, 'UniformOutput', false);

%% Sheet 1: HK_Baseline
sheetName = 'HK_Baseline';
hdr = [{'Case', 'State'}, metricNames(:)'];
rows = {};
for cIdx = 1:nCases
    for i = 1:nStates
        vals = squeeze(metricsHK(cIdx, i, :))';
        if all(isnan(vals))
            vc = repmat({'N/A'}, 1, nMetrics);
        else
            vc = num2cell(vals);
        end
        rows(end+1,:) = [{caseLabels{cIdx}, stateNames{i}}, vc]; %#ok<AGROW>
    end
end
writetable(cell2table(rows,'VariableNames',hdr), excelOut, 'Sheet', sheetName);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet 2: Per-Algorithm Comparison vs HK
for aIdx = 1:nAlgos
    algoName  = algorithmNames{aIdx};
    sheetName = sprintf('%s_vs_HK', algoName);
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    hkCols    = strcat('HK_',   metricNames);
    algoCols  = strcat(algoName,'_', metricNames);
    deltaCols = strcat('Delta_', metricNames);
    ratioCols = strcat('Ratio_', metricNames);
    hdr = [{'Case','State'}, hkCols(:)', algoCols(:)', deltaCols(:)', ratioCols(:)', {'Winner'}];

    rows = {};
    for cIdx = 1:nCases
        for i = 1:nStates
            hkVals    = squeeze(metricsHK(cIdx, i, :))';
            algoVals  = squeeze(metricsAlgo(aIdx, cIdx, i, :))';
            dVals     = squeeze(deltaMetrics(aIdx, cIdx, i, :))';
            rVals     = squeeze(ratioMetrics(aIdx, cIdx, i, :))';

            if all(isnan(algoVals))
                winner = 'N/A';
            else
                dR2 = dVals(5);
                if isnan(dR2)
                    winner = 'N/A';
                elseif dR2 > 0.001
                    winner = algoName;
                elseif dR2 < -0.001
                    winner = 'HOKALMAN';
                else
                    winner = 'Draw';
                end
            end

            toCell = @(v) num2cell(v);
            replaceNaN = @(v) mat2cell_safe(v);

            rows(end+1,:) = [{caseLabels{cIdx}, stateNames{i}}, ...
                replaceNaN(hkVals), replaceNaN(algoVals), ...
                replaceNaN(dVals),  replaceNaN(rVals), {winner}]; %#ok<AGROW>
        end
    end

    writetable(cell2table(rows,'VariableNames',hdr), excelOut, 'Sheet', sheetName);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet 3: Per-Algorithm Summary Statistics
for aIdx = 1:nAlgos
    algoName  = algorithmNames{aIdx};
    sheetName = sprintf('%s_vs_HK_Sum', algoName);
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    summHdr = {'State'};
    for m = 1:nMetrics
        summHdr{end+1} = sprintf('HK_%s_Mean',    metricNames{m}); %#ok<AGROW>
        summHdr{end+1} = sprintf('%s_%s_Mean',    algoName, metricNames{m});
        summHdr{end+1} = sprintf('Delta_%s_Mean', metricNames{m});
        summHdr{end+1} = sprintf('Delta_%s_Std',  metricNames{m});
        summHdr{end+1} = sprintf('Delta_%s_Min',  metricNames{m});
        summHdr{end+1} = sprintf('Delta_%s_Max',  metricNames{m});
    end

    summRows = {};
    for i = 1:nStates
        row = {stateNames{i}};
        for m = 1:nMetrics
            hkCol    = squeeze(metricsHK(:, i, m));
            algoCol  = squeeze(metricsAlgo(aIdx, :, i, m))';
            deltaCol = squeeze(deltaMetrics(aIdx, :, i, m))';
            row{end+1} = mean(hkCol,   'omitnan'); %#ok<AGROW>
            row{end+1} = mean(algoCol, 'omitnan');
            row{end+1} = mean(deltaCol,'omitnan');
            row{end+1} = std(deltaCol, 'omitnan');
            row{end+1} = min(deltaCol);
            row{end+1} = max(deltaCol);
        end
        summRows(end+1,:) = row; %#ok<AGROW>
    end

    row = {'[All States]'};
    for m = 1:nMetrics
        hkFlat   = reshape(metricsHK(:,:,m), [], 1);
        algoFlat = reshape(metricsAlgo(aIdx,:,:,m), [], 1);
        dFlat    = reshape(deltaMetrics(aIdx,:,:,m), [], 1);
        row{end+1} = mean(hkFlat,  'omitnan'); %#ok<AGROW>
        row{end+1} = mean(algoFlat,'omitnan');
        row{end+1} = mean(dFlat,   'omitnan');
        row{end+1} = std(dFlat,    'omitnan');
        row{end+1} = min(dFlat);
        row{end+1} = max(dFlat);
    end
    summRows(end+1,:) = row;

    writetable(cell2table(summRows,'VariableNames',summHdr), excelOut, 'Sheet', sheetName);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheets 4-9: Cross-Algorithm Comparison Matrices
crossSpecs = {
    'Cross_Delta_R2',    deltaMetrics, 5,  true;
    'Cross_Delta_FIT',   deltaMetrics, 6,  true;
    'Cross_Delta_RMSE',  deltaMetrics, 1,  false;
    'Cross_Delta_MAE',   deltaMetrics, 2,  false;
    'Cross_Delta_NRMSE', deltaMetrics, 3,  false;
    'Cross_Ratio_RMSE',  ratioMetrics, 1,  false;
    'Cross_Ratio_R2',    ratioMetrics, 5,  true;
    'Cross_Ratio_FIT',   ratioMetrics, 6,  true;
};

for rowSpec = 1:size(crossSpecs,1)
    sheetName = crossSpecs{rowSpec,1};
    dataArr   = crossSpecs{rowSpec,2};
    mIdx      = crossSpecs{rowSpec,3};
    hib       = crossSpecs{rowSpec,4};
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    valsMean = squeeze(mean(dataArr(:,:,:,mIdx), 3, 'omitnan'));

    isDelta = contains(crossSpecs{rowSpec,1},'Delta');
    isDeltaLabel = '';
    if isDelta
        if hib
            isDeltaLabel = '(+ve = algo WINS)';
        else
            isDeltaLabel = '(-ve = algo WINS)';
        end
    else
        if hib
            isDeltaLabel = '(>1 = algo WINS)';
        else
            isDeltaLabel = '(<1 = algo WINS)';
        end
    end

    hdr = [{'Algorithm'}, caseLabels(:)', ...
           {'Mean_AllCases','Std_AllCases','Best_Case','Wins_vs_HK','Losses_vs_HK'}];
    rows = {};

    for aIdx = 1:nAlgos
        row  = algorithmNames(aIdx);
        vals = valsMean(aIdx,:);
        for cIdx = 1:nCases
            if isnan(vals(cIdx)), row{end+1} = 'N/A'; %#ok<AGROW>
            else,                 row{end+1} = vals(cIdx); end
        end
        row{end+1} = mean(vals,'omitnan');
        row{end+1} = std(vals, 'omitnan');

        if hib
            [~,bi] = max(vals);
            wins   = sum(vals > 0.001,  'omitnan');
            losses = sum(vals < -0.001, 'omitnan');
        else
            [~,bi] = min(vals);
            if isDelta
                wins   = sum(vals < -0.001, 'omitnan');
                losses = sum(vals >  0.001, 'omitnan');
            else
                wins   = sum(vals < 0.999, 'omitnan');
                losses = sum(vals > 1.001, 'omitnan');
            end
        end
        row{end+1} = caseLabels{bi};
        row{end+1} = wins;
        row{end+1} = losses;
        rows(end+1,:) = row; %#ok<AGROW>
    end

    rows(end+1,:) = [{'NOTE: ' isDeltaLabel}, repmat({''}, 1, numel(hdr)-1)];

    writetable(cell2table(rows,'VariableNames',hdr), excelOut,'Sheet',sheetName);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet: Ranking_ByCase
sheetName = 'Ranking_ByCase';

maxEntries = nAlgos + 1;
hdr = {'Case'};
for r = 1:maxEntries
    hdr{end+1} = sprintf('Rank%d_Algo',  r); %#ok<AGROW>
    hdr{end+1} = sprintf('Rank%d_MeanR2',r);
    hdr{end+1} = sprintf('Rank%d_MeanFIT',r);
end
hdr{end+1} = 'HK_Rank';

rows = {};
for cIdx = 1:nCases
    cLabel = caseLabels{cIdx};

    r2vals  = zeros(1, nAlgos+1);
    fitvals = zeros(1, nAlgos+1);
    names   = [algorithmNames, {baselineName}];

    for aIdx = 1:nAlgos
        r2vals(aIdx)  = mean(squeeze(metricsAlgo(aIdx,cIdx,:,5)), 'omitnan');
        fitvals(aIdx) = mean(squeeze(metricsAlgo(aIdx,cIdx,:,6)), 'omitnan');
    end
    r2vals(end)  = mean(squeeze(metricsHK(cIdx,:,5)), 'omitnan');
    fitvals(end) = mean(squeeze(metricsHK(cIdx,:,6)), 'omitnan');

    [~, sortOrd] = sort(r2vals, 'descend');
    hkRank       = find(sortOrd == nAlgos+1);

    row = {cLabel};
    for r = 1:maxEntries
        idx = sortOrd(r);
        row{end+1} = names{idx};      %#ok<AGROW>
        row{end+1} = r2vals(idx);
        row{end+1} = fitvals(idx);
    end
    row{end+1} = hkRank;
    rows(end+1,:) = row; %#ok<AGROW>
end

writetable(cell2table(rows,'VariableNames',hdr), excelOut,'Sheet',sheetName);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet: Ranking_Overall
sheetName = 'Ranking_Overall';
hdr = {'Rank','Algorithm', ...
       'Mean_Delta_R2_AllCases',  'Std_Delta_R2', ...
       'Mean_Delta_FIT_AllCases', 'Std_Delta_FIT', ...
       'Mean_Delta_RMSE',          'Std_Delta_RMSE', ...
       'Mean_Ratio_RMSE',          'Std_Ratio_RMSE', ...
       'Mean_HK_R2',   'Mean_Algo_R2', ...
       'Mean_HK_FIT',  'Mean_Algo_FIT', ...
       'Mean_HK_RMSE', 'Mean_Algo_RMSE', ...
       'Wins_R2','Losses_R2','Draws_R2', ...
       'Overall_Verdict'};

tableData = zeros(nAlgos, 4);
for aIdx = 1:nAlgos
    dR2   = reshape(deltaMetrics(aIdx,:,:,5), [], 1);
    dFIT  = reshape(deltaMetrics(aIdx,:,:,6), [], 1);
    dRMSE = reshape(deltaMetrics(aIdx,:,:,1), [], 1);
    rRMSE = reshape(ratioMetrics(aIdx,:,:,1), [], 1);
    tableData(aIdx,:) = [mean(dR2,'omitnan'), mean(dFIT,'omitnan'), ...
                         mean(dRMSE,'omitnan'), mean(rRMSE,'omitnan')];
end

[~, sortOrd] = sort(tableData(:,1), 'descend');

rows = {};
for r = 1:nAlgos
    aIdx     = sortOrd(r);
    algoName = algorithmNames{aIdx};

    dR2   = reshape(deltaMetrics(aIdx,:,:,5), [], 1);
    dFIT  = reshape(deltaMetrics(aIdx,:,:,6), [], 1);
    dRMSE = reshape(deltaMetrics(aIdx,:,:,1), [], 1);
    rRMSE = reshape(ratioMetrics(aIdx,:,:,1), [], 1);

    hkR2flat   = reshape(metricsHK(:,:,5), [], 1);
    algoR2flat = reshape(metricsAlgo(aIdx,:,:,5), [], 1);
    hkFITflat  = reshape(metricsHK(:,:,6), [], 1);
    algoFITflat= reshape(metricsAlgo(aIdx,:,:,6), [], 1);
    hkRMSflat  = reshape(metricsHK(:,:,1), [], 1);
    algoRMSflat= reshape(metricsAlgo(aIdx,:,:,1), [], 1);

    wins=0; losses=0; draws=0;
    for cIdx = 1:nCases
        algoR2c = mean(squeeze(metricsAlgo(aIdx,cIdx,:,5)),'omitnan');
        hkR2c   = mean(squeeze(metricsHK(cIdx,:,5)),'omitnan');
        d = algoR2c - hkR2c;
        if d > 0.001,     wins   = wins+1;
        elseif d < -0.001, losses = losses+1;
        else,              draws  = draws+1;
        end
    end

    dR2Mean = mean(dR2,'omitnan');
    if dR2Mean > 0.05
        verdict = 'MUCH BETTER than HK';
    elseif dR2Mean > 0.01
        verdict = 'BETTER than HK';
    elseif dR2Mean > -0.01
        verdict = 'COMPARABLE to HK';
    elseif dR2Mean > -0.05
        verdict = 'WORSE than HK';
    else
        verdict = 'MUCH WORSE than HK';
    end

    rows(end+1,:) = {r, algoName, ...
        mean(dR2,'omitnan'),  std(dR2,'omitnan'), ...
        mean(dFIT,'omitnan'), std(dFIT,'omitnan'), ...
        mean(dRMSE,'omitnan'),std(dRMSE,'omitnan'), ...
        mean(rRMSE,'omitnan'),std(rRMSE,'omitnan'), ...
        mean(hkR2flat,'omitnan'),   mean(algoR2flat,'omitnan'), ...
        mean(hkFITflat,'omitnan'),  mean(algoFITflat,'omitnan'), ...
        mean(hkRMSflat,'omitnan'),  mean(algoRMSflat,'omitnan'), ...
        wins, losses, draws, verdict}; %#ok<AGROW>
end

writetable(cell2table(rows,'VariableNames',hdr), excelOut,'Sheet',sheetName);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet: WinsLosses
sheetName = 'WinsLosses';
hdr = [{'Algorithm'}, caseLabels(:)', ...
       {'TotalWins','TotalLosses','TotalDraws','WinRate_%'}];
rows = {};

for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    row = {algoName};
    totalW=0; totalL=0; totalD=0;
    for cIdx = 1:nCases
        algoR2c = mean(squeeze(metricsAlgo(aIdx,cIdx,:,5)),'omitnan');
        hkR2c   = mean(squeeze(metricsHK(cIdx,:,5)),'omitnan');
        d = algoR2c - hkR2c;
        if isnan(d)
            row{end+1} = 'N/A'; %#ok<AGROW>
        elseif d > 0.001
            row{end+1} = 'WIN';  totalW = totalW+1; %#ok<AGROW>
        elseif d < -0.001
            row{end+1} = 'LOSS'; totalL = totalL+1; %#ok<AGROW>
        else
            row{end+1} = 'DRAW'; totalD = totalD+1; %#ok<AGROW>
        end
    end
    validCases = totalW+totalL+totalD;
    winRate    = 100*totalW/max(validCases,1);
    row{end+1} = totalW;
    row{end+1} = totalL;
    row{end+1} = totalD;
    row{end+1} = winRate;
    rows(end+1,:) = row; %#ok<AGROW>
end

writetable(cell2table(rows,'VariableNames',hdr), excelOut,'Sheet',sheetName);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet: Per-State Delta Breakdown
for sIdx = 1:nStates
    sheetName = sprintf('Delta_%s', stateNames{sIdx});
    sheetName = regexprep(sheetName, '[()/ ]', '_');
    sheetName = regexprep(sheetName, '_+', '_');
    sheetName = regexprep(sheetName, '_$', '');
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    keyLabels = {'DeltaR2','DeltaFIT','DeltaRMSE','RatioRMSE'};
    keyIdx    = {5, 6, 1, 1};
    keyDelta  = {true, true, true, false};

    hdr = {'Algorithm'};
    for cIdx = 1:nCases
        for k = 1:numel(keyLabels)
            hdr{end+1} = sprintf('%s_%s', caseLabels{cIdx}, keyLabels{k}); %#ok<AGROW>
        end
    end
    hdr{end+1} = 'Mean_DeltaR2';
    hdr{end+1} = 'Mean_DeltaFIT';
    hdr{end+1} = 'Mean_DeltaRMSE';

    rows = {};
    for aIdx = 1:nAlgos
        row = algorithmNames(aIdx);
        for cIdx = 1:nCases
            for k = 1:numel(keyLabels)
                if keyDelta{k}
                    v = deltaMetrics(aIdx, cIdx, sIdx, keyIdx{k});
                else
                    v = ratioMetrics(aIdx, cIdx, sIdx, keyIdx{k});
                end
                if isnan(v), row{end+1} = 'N/A'; %#ok<AGROW>
                else,        row{end+1} = v; end
            end
        end
        row{end+1} = mean(squeeze(deltaMetrics(aIdx,:,sIdx,5)),'omitnan');
        row{end+1} = mean(squeeze(deltaMetrics(aIdx,:,sIdx,6)),'omitnan');
        row{end+1} = mean(squeeze(deltaMetrics(aIdx,:,sIdx,1)),'omitnan');
        rows(end+1,:) = row; %#ok<AGROW>
    end

    writetable(cell2table(rows,'VariableNames',hdr), excelOut,'Sheet',sheetName);
    fprintf('  Sheet written: %s\n', sheetName);
end

fprintf('\n========================================================\n');
fprintf('  DONE.  Output saved to:\n  %s\n', fullfile(pwd,'output',outputExcel));
fprintf('========================================================\n');
end

%% Local Functions

function [X_real, X_sim, time] = runODE(matFile)
% Simulates identified longitudinal flight dynamics model over flight test window
    S = load(matFile);

    Vb            = S.Vb;
    pqr_data      = S.pqr;
    phi_theta_psi = S.phi_theta_psi;
    Xe            = S.Xe;
    elevator      = S.elevator;
    flapPos       = S.flapPos;
    flapNeg       = S.flapNeg;
    flapDiff      = S.flapDiff;
    globalBest    = S.globalBest;
    staticData    = S.staticData;

    try, margin = S.inputData.ValidationMargin; catch, margin = 0; end

    N      = size(Vb, 1);
    dStart = 1 + margin;
    dEnd   = N - margin;

    if isfield(staticData, 'time_all')
        time = staticData.time_all;
    else
        t0   = (0:N-1)';
        time = t0(dStart:dEnd);
    end

    u_data     = Vb(dStart:dEnd, 1);
    w_data     = Vb(dStart:dEnd, 3);
    q_data     = pqr_data(dStart:dEnd, 2);
    theta_data = phi_theta_psi(2, dStart:dEnd)';
    h_data     = -Xe(dStart:dEnd, 3);
    de_data    = elevator(dStart:dEnd, 1);
    fp_data    = flapPos(dStart:dEnd, 1);
    fn_data    = flapNeg(dStart:dEnd, 1);
    fd_data    = flapDiff(dStart:dEnd, 1);

    X_real = [u_data, w_data, q_data, theta_data, h_data]';
    X0     = X_real(:, 1);
    thrust = 0;

    uo_mean = mean(u_data);
    [A_id, B_id] = formatParameters24(globalBest, uo_mean);

    odefun = @(t, x) LongODE_tv_local(t, x, A_id, B_id, ...
        time, u_data, de_data, fp_data, fn_data, fd_data, thrust);
    opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

    [t_sim, X_raw] = ode15s(odefun, time, X0, opts);
    X_sim = interp1(t_sim, X_raw, time, 'linear', 'extrap');
end

function dx = LongODE_tv_local(t, x, A, B, time, u, elevator, flapPos, flapNeg, flapDiff, thrust)
% Time-varying longitudinal state equations accounting for dynamic airspeed coupling in altitude rate
    de  = interp1(time, elevator, t, 'linear', 'extrap');
    fp  = interp1(time, flapPos,  t, 'linear', 'extrap');
    fn  = interp1(time, flapNeg,  t, 'linear', 'extrap');
    fd  = interp1(time, flapDiff, t, 'linear', 'extrap');
    ut  = interp1(time, u,        t, 'linear', 'extrap');
    U   = [de; thrust; fp; fn; fd];
    A_tv      = A;
    A_tv(5,4) = ut;
    dx  = A_tv * x + B * U;
end

function m = computeMetrics(y, yhat)
% Returns [RMSE, MAE, NRMSE, MAPE, R2, FIT] goodness-of-fit metrics
    y    = y(:);   yhat = yhat(:);
    valid = isfinite(y) & isfinite(yhat);
    y     = y(valid);  yhat = yhat(valid);
    if isempty(y), m = NaN(1,6); return; end

    err   = y - yhat;
    RMSE  = sqrt(mean(err.^2));
    MAE   = mean(abs(err));

    yRange = max(y) - min(y);
    NRMSE  = NaN;
    if yRange ~= 0, NRMSE = RMSE / yRange; end

    denom = max(abs(y), eps * max(abs(y)));
    MAPE  = 100 * mean(abs(err) ./ denom);

    SS_res = sum(err.^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = NaN; FIT = NaN;
    if SS_tot ~= 0
        R2  = 1 - SS_res / SS_tot;
        FIT = 100 * (1 - sqrt(SS_res) / sqrt(SS_tot));
    end
    m = [RMSE, MAE, NRMSE, MAPE, R2, FIT];
end

function c = mat2cell_safe(v)
% Converts numeric row vector to cell array with NaN mapped to 'N/A'
    c = cell(1, numel(v));
    for k = 1:numel(v)
        if isnan(v(k)), c{k} = 'N/A';
        else,           c{k} = v(k); end
    end
end
