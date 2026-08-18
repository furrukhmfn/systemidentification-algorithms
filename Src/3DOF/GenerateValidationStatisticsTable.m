function GenerateValidationStatisticsTable(algorithmNames, outputExcel)
% GenerateValidationStatisticsTable Computes and compiles validation statistics
% across all optimization algorithms and flight test cases for longitudinal dynamics.
%
% Validation Modes:
%   1. ODE15s State Simulation: Full time-domain integration of state equations
%      \dot{x} = A(t) x + B u over the validation window.
%   2. One-Step Derivative Prediction: Direct evaluation of model state derivatives
%      against measured kinematics (\dot{x}_{model} vs \dot{x}_{true}).
%
% Computed Metrics:
%   RMSE, MAE, NRMSE, MAPE (%), R², and FIT (%) for all longitudinal states
%   and derivative equations [u, w, q, \theta, h].
%
% Outputs:
%   Excel workbook saved to output/<outputExcel> containing detailed case sheets,
%   summary tables, cross-algorithm matrices, and overall rankings.
%
% Syntax:
%   GenerateValidationStatisticsTable()
%   GenerateValidationStatisticsTable({'ALO','GWO','CMAES'})
%   GenerateValidationStatisticsTable({'ALO','GWO'}, 'MyReport.xlsx')

%% Setup and Defaults
if nargin < 1 || isempty(algorithmNames)
    algorithmNames = {'ALO','GWO','GHOA','WOA','SSA','SCA','WCA','MFO', ...
                      'CMAES','JADE','SHADE','LSHADE'};
end
if nargin < 2 || isempty(outputExcel)
    outputExcel = 'ValidationStatisticsTable.xlsx';
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

%% State and Metric Definitions
stateNamesODE   = {'u (m/s)', 'w (m/s)', 'q (rad/s)', 'theta (rad)', 'h (m)'};
stateNames1Step = {'u_dot', 'w_dot', 'q_dot', 'theta_dot', 'h_dot'};
nStates         = 5;
metricNames     = {'RMSE', 'MAE', 'NRMSE', 'MAPE_%', 'R2', 'FIT_%'};
nMetrics        = numel(metricNames);

%% Pre-allocation
nAlgos          = numel(algorithmNames);
allMetricsODE   = NaN(nAlgos, nCases, nStates, nMetrics);
allMetrics1Step = NaN(nAlgos, nCases, nStates, nMetrics);

%% Main Processing Loop
fprintf('\n========================================================\n');
fprintf('  GenerateValidationStatisticsTable\n');
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

            try
                margin = S.inputData.ValidationMargin;
            catch
                margin = 0;
            end

            N      = size(Vb, 1);
            dStart = 1 + margin;
            dEnd   = N - margin;

            if isfield(staticData, 'time_all')
                time = staticData.time_all;
            else
                timeRaw = (0:N-1)';
                time    = timeRaw(dStart:dEnd);
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
            Nobs   = size(X_real, 2);

            uo_mean = mean(u_data);
            thrust  = 0;
            [A_id, B_id] = formatParameters24(globalBest, uo_mean);

            %% Mode 1: Full State Simulation via ode15s
            odefun = @(t, x) LongODE_tv_local(t, x, A_id, B_id, ...
                time, u_data, de_data, fp_data, fn_data, fd_data, thrust);
            opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

            [t_sim, X_raw] = ode15s(odefun, time, X0, opts);
            X_sim = interp1(t_sim, X_raw, time, 'linear', 'extrap');

            for i = 1:nStates
                y    = X_real(i, :)';
                yhat = X_sim(:, i);
                allMetricsODE(aIdx, cIdx, i, :) = computeMetrics(y, yhat);
            end

            %% Mode 2: One-Step Derivative Prediction
            U_mat       = [de_data'; thrust * ones(1, Nobs); fp_data'; fn_data'; fd_data'];
            X_dot_model = A_id * X_real + B_id * U_mat;

            if isfield(staticData, 'sim_error_check_all')
                x_dot_true = staticData.sim_error_check_all;
                if size(x_dot_true, 2) ~= Nobs
                    x_dot_true = x_dot_true(:, dStart:dEnd);
                end
            else
                dt         = mean(diff(time));
                x_dot_true = zeros(5, Nobs);
                x_dot_true(:, 2:end-1) = (X_real(:,3:end) - X_real(:,1:end-2)) / (2*dt);
                x_dot_true(:, 1)       = (X_real(:,2) - X_real(:,1)) / dt;
                x_dot_true(:, end)     = (X_real(:,end) - X_real(:,end-1)) / dt;
            end

            for i = 1:nStates
                y    = x_dot_true(i, :)';
                yhat = X_dot_model(i, :)';
                allMetrics1Step(aIdx, cIdx, i, :) = computeMetrics(y, yhat);
            end

            fprintf('  Case %d (%d/%d) - OK\n', caseNum, cIdx, nCases);

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

%% Sheet A: Per-Algorithm Validation Detail
for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};

    for modeIdx = 1:2
        if modeIdx == 1
            metrics   = squeeze(allMetricsODE(aIdx, :, :, :));
            sheetName = sprintf('%s_ODE_Detail', algoName);
            sNames    = stateNamesODE;
        else
            metrics   = squeeze(allMetrics1Step(aIdx, :, :, :));
            sheetName = sprintf('%s_1Step_Detail', algoName);
            sNames    = stateNames1Step;
        end
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
                dataRows(end+1, :) = [{cLabel, sNames{i}}, valCell]; %#ok<AGROW>
            end
        end

        T = cell2table(dataRows, 'VariableNames', headerRow);
        writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
        fprintf('  Sheet written: %s\n', sheetName);
    end
end

%% Sheet B: Per-Algorithm Summary Statistics
for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};

    for modeIdx = 1:2
        if modeIdx == 1
            metrics   = squeeze(allMetricsODE(aIdx, :, :, :));
            sheetName = sprintf('%s_ODE_Summary', algoName);
            sNames    = stateNamesODE;
        else
            metrics   = squeeze(allMetrics1Step(aIdx, :, :, :));
            sheetName = sprintf('%s_1Step_Summary', algoName);
            sNames    = stateNames1Step;
        end
        if numel(sheetName) > 31, sheetName = sheetName(1:31); end

        summaryHeader = {'State'};
        for m = 1:nMetrics
            summaryHeader{end+1} = [metricNames{m} '_Mean']; %#ok<AGROW>
            summaryHeader{end+1} = [metricNames{m} '_Std'];
            summaryHeader{end+1} = [metricNames{m} '_Min'];
            summaryHeader{end+1} = [metricNames{m} '_Max'];
        end

        summaryRows = {};
        for i = 1:nStates
            vals_all = squeeze(metrics(:, i, :));
            row = {sNames{i}};
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

        T = cell2table(summaryRows, 'VariableNames', summaryHeader);
        writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
        fprintf('  Sheet written: %s\n', sheetName);
    end
end

%% Sheet C: Cross-Algorithm Metric Comparisons
crossMetrics    = {'R2', 'FIT_pct', 'RMSE', 'MAE', 'NRMSE', 'MAPE_pct'};
crossMetricIdx  = [5, 6, 1, 2, 3, 4];
higherIsBetter  = [true, true, false, false, false, false];

for modeIdx = 1:2
    modeLbl = {'ODE', '1Step'};
    if modeIdx == 1
        allM = allMetricsODE;
    else
        allM = allMetrics1Step;
    end

    for cmIdx = 1:numel(crossMetrics)
        mName     = crossMetrics{cmIdx};
        mIdx      = crossMetricIdx(cmIdx);
        sheetName = sprintf('Cross_%s_%s', modeLbl{modeIdx}, mName);
        if numel(sheetName) > 31, sheetName = sheetName(1:31); end

        valsMean  = squeeze(mean(allM(:, :, :, mIdx), 3, 'omitnan'));

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
end

%% Sheet D: Best Algorithm per Case
sheetName = 'BestAlgo_PerCase';
headerRow = {'Case', 'Best_Algo_ODE_R2', 'Best_R2_ODE', ...
             'Best_Algo_1Step_R2', 'Best_R2_1Step', ...
             'Best_Algo_ODE_FIT',  'Best_FIT_ODE'};
dataRows  = {};

for cIdx = 1:nCases
    cLabel   = caseLabels{cIdx};
    r2ode    = squeeze(mean(allMetricsODE(:,  cIdx, :, 5), 3, 'omitnan'));
    r21s     = squeeze(mean(allMetrics1Step(:, cIdx, :, 5), 3, 'omitnan'));
    fitode   = squeeze(mean(allMetricsODE(:,  cIdx, :, 6), 3, 'omitnan'));

    [bestR2ode,  iOde]  = max(r2ode);
    [bestR21s,   i1s]   = max(r21s);
    [bestFITode, iFit]  = max(fitode);

    dataRows(end+1, :) = {cLabel, ...
        algorithmNames{iOde},  bestR2ode,  ...
        algorithmNames{i1s},   bestR21s,   ...
        algorithmNames{iFit},  bestFITode}; %#ok<AGROW>
end

T = cell2table(dataRows, 'VariableNames', headerRow);
writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet E: Grand Validation Summary
sheetName  = 'GrandSummary';
nM         = nMetrics;
hdrMean    = strcat(metricNames, '_Mean');
hdrStd     = strcat(metricNames, '_Std');
sumHdr     = [{'Algorithm', 'Mode'}, hdrMean(:)', hdrStd(:)'];

grandRows  = {};
for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    for modeIdx = 1:2
        if modeIdx == 1
            allM    = allMetricsODE;
            modeLbl = 'ODE';
        else
            allM    = allMetrics1Step;
            modeLbl = '1Step';
        end
        flat     = reshape(squeeze(allM(aIdx, :, :, :)), nCases*nStates, nM);
        meanVals = mean(flat, 1, 'omitnan');
        stdVals  = std(flat,  0, 1, 'omitnan');
        grandRows(end+1, :) = [{algoName, modeLbl}, num2cell(meanVals), num2cell(stdVals)]; %#ok<AGROW>
    end
end

T = cell2table(grandRows, 'VariableNames', sumHdr);
writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
fprintf('  Sheet written: %s\n', sheetName);

%% Sheet F: Per-State Detailed Comparisons
for sIdx = 1:nStates
    for modeIdx = 1:2
        if modeIdx == 1
            allM      = allMetricsODE;
            sLabel    = stateNamesODE{sIdx};
            modeLbl   = 'ODE';
        else
            allM      = allMetrics1Step;
            sLabel    = stateNames1Step{sIdx};
            modeLbl   = '1Step';
        end

        sheetName = sprintf('%s_%s', modeLbl, sLabel);
        sheetName = regexprep(sheetName, '[()/ ]', '_');
        sheetName = regexprep(sheetName, '_+', '_');
        sheetName = regexprep(sheetName, '_$', '');
        if numel(sheetName) > 31, sheetName = sheetName(1:31); end

        keyMetrics   = {'R2', 'FIT_%', 'RMSE', 'MAE'};
        keyMetricIdx = [5, 6, 1, 2];

        headerCols = {'Algorithm'};
        for cIdx = 1:nCases
            for km = 1:numel(keyMetrics)
                headerCols{end+1} = sprintf('%s_%s', caseLabels{cIdx}, keyMetrics{km}); %#ok<AGROW>
            end
        end
        headerCols{end+1} = 'Mean_R2';
        headerCols{end+1} = 'Mean_FIT';
        headerCols{end+1} = 'Mean_RMSE';

        dataRows = {};
        for aIdx = 1:nAlgos
            row = algorithmNames(aIdx);
            for cIdx = 1:nCases
                for km = 1:numel(keyMetrics)
                    v = allM(aIdx, cIdx, sIdx, keyMetricIdx(km));
                    if isnan(v)
                        row{end+1} = 'N/A'; %#ok<AGROW>
                    else
                        row{end+1} = v; %#ok<AGROW>
                    end
                end
            end
            row{end+1} = mean(squeeze(allM(aIdx, :, sIdx, 5)), 'omitnan');
            row{end+1} = mean(squeeze(allM(aIdx, :, sIdx, 6)), 'omitnan');
            row{end+1} = mean(squeeze(allM(aIdx, :, sIdx, 1)), 'omitnan');
            dataRows(end+1, :) = row; %#ok<AGROW>
        end

        T = cell2table(dataRows, 'VariableNames', headerCols);
        writetable(T, excelOut, 'Sheet', sheetName, 'WriteVariableNames', true);
        fprintf('  Sheet written: %s\n', sheetName);
    end
end

fprintf('\n========================================================\n');
fprintf('  DONE.  Output saved to:\n  %s\n', fullfile(pwd, 'output', outputExcel));
fprintf('========================================================\n');
end

%% Local Helper Functions

function m = computeMetrics(y, yhat)
% Returns [RMSE, MAE, NRMSE, MAPE, R2, FIT] goodness-of-fit metrics
    y    = y(:);
    yhat = yhat(:);
    valid = isfinite(y) & isfinite(yhat);
    y     = y(valid);
    yhat  = yhat(valid);

    if isempty(y)
        m = NaN(1, 6);
        return
    end

    err  = y - yhat;

    RMSE  = sqrt(mean(err.^2));
    MAE   = mean(abs(err));

    yRange = max(y) - min(y);
    if yRange == 0
        NRMSE = NaN;
    else
        NRMSE = RMSE / yRange;
    end

    denom = max(abs(y), eps * max(abs(y)));
    MAPE  = 100 * mean(abs(err) ./ denom);

    SS_res = sum(err.^2);
    SS_tot = sum((y - mean(y)).^2);
    if SS_tot == 0
        R2  = NaN;
        FIT = NaN;
    else
        R2  = 1 - SS_res / SS_tot;
        FIT = 100 * (1 - sqrt(SS_res) / sqrt(SS_tot));
    end

    m = [RMSE, MAE, NRMSE, MAPE, R2, FIT];
end

function dx = LongODE_tv_local(t, x, A, B, time, u, elevator, flapPos, flapNeg, flapDiff, thrust)
% Longitudinal state derivative with time-varying airspeed u(t) for altitude rate
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
