function GenerateTrainingStatisticsTable(algorithmNames, outputExcel)
% GenerateTrainingStatisticsTable Compiles a comprehensive statistical summary
% of optimizer training performance, convergence characteristics, and derivative
% fitting quality across all metaheuristic algorithms and 20 flight test cases.
%
% Extracted Metrics:
%   - Convergence: Initial/final RMSE, percentage improvement, completed iterations,
%     best-update iteration index, population statistics, and convergence slope.
%   - Per-State Derivative Fit: RMSE, MAE, NRMSE, MAPE, R², and FIT% for longitudinal
%     derivatives (u_dot, w_dot, q_dot, theta_dot, h_dot).
%   - Weighted Objective: Final best-fitness value and state weight vectors.
%
% Outputs:
%   Excel report written to output/<outputExcel> containing per-algorithm details,
%   cross-algorithm comparisons, case rankings, and iteration learning curves.
%
% Syntax:
%   GenerateTrainingStatisticsTable()
%   GenerateTrainingStatisticsTable({'ALO','GWO','CMAES'})
%   GenerateTrainingStatisticsTable({'ALO','GWO'}, 'MyTrainingReport.xlsx')

%% Setup and Defaults
if nargin < 1 || isempty(algorithmNames)
    algorithmNames = {'ALO','GWO','GHOA','WOA','SSA','SCA','WCA','MFO', ...
                      'CMAES','JADE','SHADE','LSHADE'};
end
if nargin < 2 || isempty(outputExcel)
    outputExcel = 'TrainingStatisticsTable.xlsx';
end

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

%% Field and Metric Definitions
stateNames   = {'u_dot', 'w_dot', 'q_dot', 'theta_dot', 'h_dot'};
nStates      = 5;
metricNames  = {'RMSE', 'MAE', 'NRMSE', 'MAPE_%', 'R2', 'FIT_%'};
nMetrics     = numel(metricNames);

convNames = {'InitialRMSE', 'FinalWeightedRMSE', 'ImprovementPct', ...
             'IterationsRan', 'IterBestFound', ...
             'FinalMeanRMSE', 'FinalStdRMSE', 'FinalMinRMSE', 'FinalMaxRMSE', ...
             'ConvergenceSlope', 'PopulationSize', 'MaxIterations', ...
             'W_u', 'W_w', 'W_q', 'W_theta', 'W_h'};
nConv        = numel(convNames);

%% Pre-allocation
nAlgos          = numel(algorithmNames);
allMetricsTrn   = NaN(nAlgos, nCases, nStates, nMetrics);
allConvStats    = NaN(nAlgos, nCases, nConv);

%% Main Processing Loop
fprintf('\n========================================================\n');
fprintf('  GenerateTrainingStatisticsTable\n');
fprintf('  Algorithms : %d   |  Cases per algo : %d\n', nAlgos, nCases);
fprintf('========================================================\n');

for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    fprintf('\n[%d/%d] Processing algorithm: %s\n', aIdx, nAlgos, algoName);

    for cIdx = 1:nCases
        caseNum = caseNums(cIdx);
        matFile = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, caseNum));

        if ~exist(matFile, 'file')
            fprintf('  [SKIP] %s not found.\n', matFile);
            continue
        end

        try
            S = load(matFile);

            globalBest        = S.globalBest;
            globalBestFitness = S.globalBestFitness;
            globalFitnessMatrix = S.globalFitnessMatrix;
            staticData        = S.staticData;
            iter              = S.iter;
            populationSize    = S.populationSize;
            maxIterations     = S.maxIterations;

            if isfield(S, 'globalrmse')
                globalrmse = S.globalrmse;
                globalrmse = globalrmse(1:iter);
            else
                globalrmse = mean(globalFitnessMatrix(1:iter, :), 2);
            end

            fitnessMatrix = globalFitnessMatrix(1:iter, :);

            if isfield(staticData, 'Weights')
                W = staticData.Weights(:)';
            else
                W = [0.5, 0.2, 0.2, 0.05, 0.05];
            end

            %% Section 1: Optimizer Convergence Statistics
            InitialRMSE     = mean(fitnessMatrix(1, :));
            FinalWeightedRMSE = globalBestFitness;
            ImprovementPct  = (InitialRMSE - FinalWeightedRMSE) / InitialRMSE * 100;
            IterationsRan   = iter;

            bestPerIter     = min(fitnessMatrix, [], 2);
            [~, IterBestFound] = min(bestPerIter);

            FinalMeanRMSE   = mean(fitnessMatrix(end, :));
            FinalStdRMSE    = std(fitnessMatrix(end, :));
            FinalMinRMSE    = min(fitnessMatrix(end, :));
            FinalMaxRMSE    = max(fitnessMatrix(end, :));

            itVec = (1:iter)';
            p     = polyfit(itVec, globalrmse(:), 1);
            ConvergenceSlope = p(1);

            allConvStats(aIdx, cIdx, :) = [InitialRMSE, FinalWeightedRMSE, ImprovementPct, ...
                IterationsRan, IterBestFound, ...
                FinalMeanRMSE, FinalStdRMSE, FinalMinRMSE, FinalMaxRMSE, ...
                ConvergenceSlope, populationSize, maxIterations, W];

            %% Section 2: Per-State Derivative Fit
            u_all     = staticData.u_all;
            w_all     = staticData.w_all;
            q_all     = staticData.q_all;
            theta_all = staticData.theta_all;
            elev_all  = staticData.elevator_all;
            fp_all    = staticData.flapPos_all;
            fn_all    = staticData.flapNeg_all;
            fd_all    = staticData.flapDiff_all;
            thrust    = staticData.thrust;
            g         = staticData.valueOfGravitationConstant;

            x_dot_true = staticData.sim_error_check_all;
            N          = size(u_all, 1);

            p1  = globalBest(1);  p2  = globalBest(2);  p3  = globalBest(3);
            p4  = globalBest(4);  p5  = globalBest(5);  p6  = globalBest(6);
            p7  = globalBest(7);  p8  = globalBest(8);  p9  = globalBest(9);
            p10 = globalBest(10); p11 = globalBest(11);
            p12 = globalBest(12); p13 = globalBest(13);
            p14 = globalBest(14); p15 = globalBest(15);
            p16 = globalBest(16); p17 = globalBest(17); p18 = globalBest(18);
            p19 = globalBest(19); p20 = globalBest(20); p21 = globalBest(21);
            p22 = globalBest(22); p23 = globalBest(23); p24 = globalBest(24);

            x_dot_pred = zeros(5, N);
            x_dot_pred(1,:) = p1*u_all' + p2*w_all' + p3*q_all' - g*theta_all' ...
                            + p10*elev_all' + p11*thrust + p16*fp_all' + p17*fn_all' + p18*fd_all';
            x_dot_pred(2,:) = p4*u_all' + p5*w_all' + p6*q_all' ...
                            + p12*elev_all' + p13*thrust + p19*fp_all' + p20*fn_all' + p21*fd_all';
            x_dot_pred(3,:) = p7*u_all' + p8*w_all' + p9*q_all' ...
                            + p14*elev_all' + p15*thrust + p22*fp_all' + p23*fn_all' + p24*fd_all';
            x_dot_pred(4,:) = q_all';
            x_dot_pred(5,:) = -w_all' + u_all' .* theta_all';

            for i = 1:nStates
                y    = x_dot_true(i, :)';
                yhat = x_dot_pred(i, :)';
                allMetricsTrn(aIdx, cIdx, i, :) = computeMetrics(y, yhat);
            end

            fprintf('  Case %d (%d/%d) - OK  [iter=%d, bestFit=%.4f]\n', ...
                caseNum, cIdx, nCases, iter, globalBestFitness);

        catch ME
            fprintf('  Case %d (%d/%d) - ERROR: %s\n', caseNum, cIdx, nCases, ME.message);
        end
    end
end

%% Excel Export Setup
if ~exist('output', 'dir'), mkdir('output'); end
excelOut = fullfile('output', outputExcel);
if exist(excelOut, 'file'), delete(excelOut); end

fprintf('\n\nWriting results to: %s\n', excelOut);

caseLabels = arrayfun(@(n) sprintf('(%d,%d)', floor(n/10), mod(n,10)), ...
    caseNums, 'UniformOutput', false);

%% Sheet A: Per-Algorithm Training Detail
for aIdx = 1:nAlgos
    algoName  = algorithmNames{aIdx};
    metrics   = squeeze(allMetricsTrn(aIdx, :, :, :));
    sheetName = sprintf('%s_Trn_Detail', algoName);
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    headerRow = [{'Case', 'State'}, metricNames(:)'];
    dataRows  = {};

    for cIdx = 1:nCases
        cLabel = caseLabels{cIdx};
        for i = 1:nStates
            vals = squeeze(metrics(cIdx, i, :))';
            if all(isnan(vals))
                valCell = repmat({'N/A'}, 1, nMetrics);
            else
                valCell = num2cell(vals);
            end
            dataRows(end+1, :) = [{cLabel, stateNames{i}}, valCell]; %#ok<AGROW>
        end
    end

    T = cell2table(dataRows, 'VariableNames', headerRow);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet B: Per-Algorithm Training Summary
for aIdx = 1:nAlgos
    algoName  = algorithmNames{aIdx};
    metrics   = squeeze(allMetricsTrn(aIdx, :, :, :));
    sheetName = sprintf('%s_Trn_Summary', algoName);
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    summaryHdr = {'State'};
    for m = 1:nMetrics
        summaryHdr{end+1} = [metricNames{m} '_Mean']; %#ok<AGROW>
        summaryHdr{end+1} = [metricNames{m} '_Std'];
        summaryHdr{end+1} = [metricNames{m} '_Min'];
        summaryHdr{end+1} = [metricNames{m} '_Max'];
    end

    summaryRows = {};
    for i = 1:nStates
        vals_all = squeeze(metrics(:, i, :));
        row = {stateNames{i}};
        for m = 1:nMetrics
            col = vals_all(:, m);
            row{end+1} = mean(col, 'omitnan'); %#ok<AGROW>
            row{end+1} = std(col,  'omitnan');
            row{end+1} = min(col);
            row{end+1} = max(col);
        end
        summaryRows(end+1, :) = row; %#ok<AGROW>
    end

    allVals = reshape(metrics, nCases*nStates, nMetrics);
    row = {'[All States]'};
    for m = 1:nMetrics
        col = allVals(:, m);
        row{end+1} = mean(col, 'omitnan'); %#ok<AGROW>
        row{end+1} = std(col,  'omitnan');
        row{end+1} = min(col);
        row{end+1} = max(col);
    end
    summaryRows(end+1, :) = row;

    T = cell2table(summaryRows, 'VariableNames', summaryHdr);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet C: Per-Algorithm Convergence Detail
for aIdx = 1:nAlgos
    algoName  = algorithmNames{aIdx};
    sheetName = sprintf('%s_Conv_Detail', algoName);
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    headerRow = [{'Case'}, convNames(:)'];
    dataRows  = {};

    for cIdx = 1:nCases
        cLabel = caseLabels{cIdx};
        vals   = squeeze(allConvStats(aIdx, cIdx, :))';
        if all(isnan(vals))
            valCell = repmat({'N/A'}, 1, nConv);
        else
            valCell = num2cell(vals);
        end
        dataRows(end+1, :) = [{cLabel}, valCell]; %#ok<AGROW>
    end

    convData = squeeze(allConvStats(aIdx, :, :));
    meanRow  = {'[Mean]'};
    stdRow   = {'[Std]'};
    minRow   = {'[Min]'};
    maxRow   = {'[Max]'};
    for k = 1:nConv
        col      = convData(:, k);
        meanRow{end+1} = mean(col, 'omitnan'); %#ok<AGROW>
        stdRow{end+1}  = std(col,  'omitnan');
        minRow{end+1}  = min(col);
        maxRow{end+1}  = max(col);
    end
    dataRows(end+1, :) = meanRow;
    dataRows(end+1, :) = stdRow;
    dataRows(end+1, :) = minRow;
    dataRows(end+1, :) = maxRow;

    T = cell2table(dataRows, 'VariableNames', headerRow);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet D: Cross-Algorithm Metric Comparisons
crossMetrics    = {'R2', 'FIT_pct', 'RMSE', 'MAE', 'NRMSE', 'MAPE_pct'};
crossMetricIdx  = [5, 6, 1, 2, 3, 4];
higherIsBetter  = [true, true, false, false, false, false];

for cmIdx = 1:numel(crossMetrics)
    mName    = crossMetrics{cmIdx};
    mIdx     = crossMetricIdx(cmIdx);
    sheetName = sprintf('Cross_Trn_%s', mName);
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    valsMean = squeeze(mean(allMetricsTrn(:, :, :, mIdx), 3, 'omitnan'));

    headerRow = [{'Algorithm'}, caseLabels(:)', {'Mean_AllCases', 'Std_AllCases', 'Best_Case'}];
    dataRows  = {};

    for aIdx = 1:nAlgos
        row  = algorithmNames(aIdx);
        vals = valsMean(aIdx, :);
        for cIdx = 1:nCases
            if isnan(vals(cIdx))
                row{end+1} = 'N/A'; %#ok<AGROW>
            else
                row{end+1} = vals(cIdx); %#ok<AGROW>
            end
        end
        row{end+1} = mean(vals, 'omitnan');
        row{end+1} = std(vals,  'omitnan');
        if higherIsBetter(cmIdx)
            [~, bestCI] = max(vals);
        else
            [~, bestCI] = min(vals);
        end
        row{end+1} = caseLabels{bestCI};
        dataRows(end+1, :) = row; %#ok<AGROW>
    end

    T = cell2table(dataRows, 'VariableNames', headerRow);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet E: Cross-Algorithm Convergence Comparisons
convCross = {
    'Cross_Conv_FinalRMSE',    3,  false;
    'Cross_Conv_ImprovPct',    4,  true;
    'Cross_Conv_Iterations',   5,  false;
    'Cross_Conv_InitRMSE',     2,  false;
    'Cross_Conv_IterBest',     6,  false;
    'Cross_Conv_FinalMean',    7,  false;
    'Cross_Conv_ConvSlope',    11, false;
};

for rowIdx = 1:size(convCross, 1)
    sheetName = convCross{rowIdx, 1};
    statIdx   = convCross{rowIdx, 2};
    hib       = convCross{rowIdx, 3};
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    valsMat = squeeze(allConvStats(:, :, statIdx));

    headerRow = [{'Algorithm'}, caseLabels(:)', {'Mean_AllCases', 'Std_AllCases', 'Best_Case'}];
    dataRows  = {};

    for aIdx = 1:nAlgos
        row  = algorithmNames(aIdx);
        vals = valsMat(aIdx, :);
        for cIdx = 1:nCases
            if isnan(vals(cIdx))
                row{end+1} = 'N/A'; %#ok<AGROW>
            else
                row{end+1} = vals(cIdx); %#ok<AGROW>
            end
        end
        row{end+1} = mean(vals, 'omitnan');
        row{end+1} = std(vals,  'omitnan');
        if hib
            [~, bestCI] = max(vals);
        else
            [~, bestCI] = min(vals);
        end
        row{end+1} = caseLabels{bestCI};
        dataRows(end+1, :) = row; %#ok<AGROW>
    end

    T = cell2table(dataRows, 'VariableNames', headerRow);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet F: Best-Performing Algorithm per Case
sheetName = 'BestAlgo_Trn_PerCase';
headerRow = {'Case', ...
             'BestAlgo_R2_Train',   'Best_R2_Train', ...
             'BestAlgo_FIT_Train',  'Best_FIT_Train', ...
             'BestAlgo_RMSE_Train', 'Best_RMSE_Train', ...
             'BestAlgo_WeightedRMSE', 'Best_WeightedRMSE'};
dataRows  = {};

for cIdx = 1:nCases
    cLabel = caseLabels{cIdx};

    r2vals   = squeeze(mean(allMetricsTrn(:, cIdx, :, 5), 3, 'omitnan'));
    fitvals  = squeeze(mean(allMetricsTrn(:, cIdx, :, 6), 3, 'omitnan'));
    rmsvals  = squeeze(mean(allMetricsTrn(:, cIdx, :, 1), 3, 'omitnan'));
    wrmse    = squeeze(allConvStats(:, cIdx, 3));

    [bestR2,   iR2]  = max(r2vals);
    [bestFIT,  iFIT] = max(fitvals);
    [bestRMSE, iRMS] = min(rmsvals);
    [bestWRMSE,iWRM] = min(wrmse);

    dataRows(end+1, :) = {cLabel, ...
        algorithmNames{iR2},  bestR2, ...
        algorithmNames{iFIT}, bestFIT, ...
        algorithmNames{iRMS}, bestRMSE, ...
        algorithmNames{iWRM}, bestWRMSE}; %#ok<AGROW>
end

T = cell2table(dataRows, 'VariableNames', headerRow);
writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet G: Grand Training Summary
sheetName = 'GrandTrainingSummary';

hdrFit   = strcat(metricNames, '_Mean');
hdrFitS  = strcat(metricNames, '_Std');
convSumCols = {'FinalWeightedRMSE_Mean', 'FinalWeightedRMSE_Std', ...
               'ImprovePct_Mean',         'ImprovePct_Std', ...
               'IterationsRan_Mean',       'IterationsRan_Std', ...
               'IterBestFound_Mean',        'IterBestFound_Std', ...
               'ConvSlope_Mean',            'ConvSlope_Std'};

grandHdr  = [{'Algorithm'}, hdrFit(:)', hdrFitS(:)', convSumCols];
grandRows = {};

convIdxMap = [2, 3, 4, 5, 10];

for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    flat     = reshape(squeeze(allMetricsTrn(aIdx, :, :, :)), nCases*nStates, nMetrics);
    meanFit  = mean(flat, 1, 'omitnan');
    stdFit   = std(flat,  0, 1, 'omitnan');

    convFlat = squeeze(allConvStats(aIdx, :, convIdxMap));
    convVals = [];
    for k = 1:numel(convIdxMap)
        col = convFlat(:, k);
        convVals = [convVals, mean(col,'omitnan'), std(col,'omitnan')]; %#ok<AGROW>
    end

    grandRows(end+1, :) = [{algoName}, num2cell(meanFit), num2cell(stdFit), num2cell(convVals)]; %#ok<AGROW>
end

fitColInGrand = 6;
if ~isempty(grandRows)
    fitScores = cellfun(@(r) r{1 + fitColInGrand}, grandRows, 'UniformOutput', false);
    fitScores = cell2mat(fitScores);
    [~, sortOrd] = sort(fitScores, 'descend');
    grandRows   = grandRows(sortOrd, :);
end

T = cell2table(grandRows, 'VariableNames', grandHdr);
writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet H: Per-State Training Detail
for sIdx = 1:nStates
    sheetName = sprintf('Trn_%s', stateNames{sIdx});
    if numel(sheetName) > 31, sheetName = sheetName(1:31); end

    keyMetrics    = {'R2', 'FIT_%', 'RMSE', 'MAE'};
    keyMetricIdx  = [5, 6, 1, 2];

    headerCols = {'Algorithm'};
    for cIdx = 1:nCases
        for km = 1:numel(keyMetrics)
            headerCols{end+1} = sprintf('%s_%s', caseLabels{cIdx}, keyMetrics{km}); %#ok<AGROW>
        end
    end
    headerCols{end+1} = 'Mean_R2';
    headerCols{end+1} = 'Mean_FIT';
    headerCols{end+1} = 'Mean_RMSE';
    headerCols{end+1} = 'Mean_MAE';

    dataRows = {};
    for aIdx = 1:nAlgos
        row = algorithmNames(aIdx);
        for cIdx = 1:nCases
            for km = 1:numel(keyMetrics)
                v = allMetricsTrn(aIdx, cIdx, sIdx, keyMetricIdx(km));
                if isnan(v)
                    row{end+1} = 'N/A'; %#ok<AGROW>
                else
                    row{end+1} = v; %#ok<AGROW>
                end
            end
        end
        row{end+1} = mean(squeeze(allMetricsTrn(aIdx, :, sIdx, 5)), 'omitnan');
        row{end+1} = mean(squeeze(allMetricsTrn(aIdx, :, sIdx, 6)), 'omitnan');
        row{end+1} = mean(squeeze(allMetricsTrn(aIdx, :, sIdx, 1)), 'omitnan');
        row{end+1} = mean(squeeze(allMetricsTrn(aIdx, :, sIdx, 2)), 'omitnan');
        dataRows(end+1, :) = row; %#ok<AGROW>
    end

    T = cell2table(dataRows, 'VariableNames', headerCols);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

%% Sheet I: Convergence Learning Curves
sheetName = 'Conv_LearningCurves';
if numel(sheetName) > 31, sheetName = sheetName(1:31); end

maxIter = 0;
for aIdx = 1:nAlgos
    for cIdx = 1:nCases
        v = allConvStats(aIdx, cIdx, 5);
        if ~isnan(v) && v > maxIter
            maxIter = v;
        end
    end
end

if maxIter > 0
    headerCols = {'Iteration'};
    for aIdx = 1:nAlgos
        for cIdx = 1:nCases
            headerCols{end+1} = sprintf('%s_%s_MeanRMSE', algorithmNames{aIdx}, caseLabels{cIdx}); %#ok<AGROW>
        end
    end

    curveData = NaN(maxIter, nAlgos * nCases);
    colIdx    = 0;
    for aIdx = 1:nAlgos
        for cIdx = 1:nCases
            colIdx  = colIdx + 1;
            caseNum = caseNums(cIdx);
            matFile = fullfile('output', sprintf('%s_Case%d_24.mat', algorithmNames{aIdx}, caseNum));
            if ~exist(matFile, 'file'), continue; end
            try
                S = load(matFile, 'globalrmse', 'globalFitnessMatrix', 'iter');
                actualIter = S.iter;
                if isfield(S, 'globalrmse') && ~isempty(S.globalrmse)
                    curve = S.globalrmse(1:actualIter);
                else
                    curve = mean(S.globalFitnessMatrix(1:actualIter, :), 2);
                end
                curveData(1:actualIter, colIdx) = curve(:);
            catch
            end
        end
    end

    iterCol  = (1:maxIter)';
    fullData = [iterCol, curveData];

    dataRows = {};
    for r = 1:maxIter
        dataRows(end+1, :) = num2cell(fullData(r, :)); %#ok<AGROW>
    end
    T = cell2table(dataRows, 'VariableNames', headerCols);
    writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
    fprintf('  Sheet written: %s\n', sheetName);
end

fprintf('\n========================================================\n');
fprintf('  DONE.  Output saved to:\n  %s\n', fullfile(pwd, 'output', outputExcel));
fprintf('========================================================\n');
end

%% Local Functions

function m = computeMetrics(y, yhat)
% Returns [RMSE, MAE, NRMSE, MAPE, R2, FIT] for flight state vectors
    y    = y(:);
    yhat = yhat(:);
    valid = isfinite(y) & isfinite(yhat);
    y     = y(valid);
    yhat  = yhat(valid);

    if isempty(y)
        m = NaN(1, 6);
        return
    end

    err   = y - yhat;
    RMSE  = sqrt(mean(err.^2));
    MAE   = mean(abs(err));

    yRange = max(y) - min(y);
    NRMSE  = NaN;
    if yRange ~= 0
        NRMSE = RMSE / yRange;
    end

    denom = max(abs(y), eps * max(abs(y)));
    MAPE  = 100 * mean(abs(err) ./ denom);

    SS_res = sum(err.^2);
    SS_tot = sum((y - mean(y)).^2);
    R2  = NaN;
    FIT = NaN;
    if SS_tot ~= 0
        R2  = 1 - SS_res / SS_tot;
        FIT = 100 * (1 - sqrt(SS_res) / sqrt(SS_tot));
    end

    m = [RMSE, MAE, NRMSE, MAPE, R2, FIT];
end
