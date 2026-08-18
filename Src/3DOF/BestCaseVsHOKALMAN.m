function BestCaseVsHOKALMAN(algorithmNames, outputExcel)
% BestCaseVsHOKALMAN Identifies the highest-performing case for each metaheuristic
% optimizer (based on mean R²) and compares it against the Ho-Kalman baseline
% across 5 longitudinal flight dynamics evaluation metrics:
%
%   1. Accuracy         - Mean Absolute Error (MAE) across all 5 states.
%   2. RMSE             - Root Mean Squared Error across all 5 states.
%   3. Stability        - Standard deviation of per-state RMSE across the test cases.
%   4. Convergence Speed- Iteration fraction (0 to 1) required to reach within 1%
%                         of the final optimal fitness.
%   5. R²               - Coefficient of Determination across all 5 states.
%
% Results are printed to the console and exported to an Excel workbook.
%
% Syntax:
%   BestCaseVsHOKALMAN()
%   BestCaseVsHOKALMAN({'ALO','GWO','CMAES'})
%   BestCaseVsHOKALMAN({'ALO','GWO'}, 'MyReport.xlsx')

%% Setup and Defaults
if nargin < 1 || isempty(algorithmNames)
    algorithmNames = {'ALO','GWO','GHOA','WOA','SSA','SCA','WCA','MFO', ...
                      'CMAES','JADE','SHADE','LSHADE'};
end
if nargin < 2 || isempty(outputExcel)
    outputExcel = 'BestCase_vs_HOKALMAN.xlsx';
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
caseLabels = arrayfun(@(n) sprintf('(%d,%d)', floor(n/10), mod(n,10)), ...
    caseNums, 'UniformOutput', false);

nAlgos = numel(algorithmNames);

%% Path Configuration
algoDir = fullfile(fileparts(mfilename('fullpath')), 'more_algorithms');
if exist(algoDir, 'dir'), addpath(algoDir); end

%% Execution Header
fprintf('\n================================================================\n');
fprintf('  BestCaseVsHOKALMAN\n');
fprintf('  Algorithms : %d    Cases : %d each\n', nAlgos, nCases);
fprintf('================================================================\n');

%% Baseline Ho-Kalman Evaluation
fprintf('\n--- Step 1: Computing Ho-Kalman baseline ---\n');
hkRMSE  = NaN(nCases, 5);   % [nCases x 5 states]
hkMAE   = NaN(nCases, 5);
hkR2    = NaN(nCases, 5);
hkConverge = NaN(nCases, 1);   % algebraic HK method has instant convergence (0)

for cIdx = 1:nCases
    caseNum = caseNums(cIdx);
    hkFile  = fullfile('output', sprintf('HOKALMAN_Case%d_24.mat', caseNum));
    if ~exist(hkFile, 'file')
        fprintf('  [SKIP HK] %s not found.\n', hkFile);
        continue
    end
    try
        [X_real, X_sim, ~] = runODE(hkFile);
        for s = 1:5
            [rmse, mae, r2] = metrics3(X_real(s,:)', X_sim(:,s));
            hkRMSE(cIdx, s) = rmse;
            hkMAE(cIdx,  s) = mae;
            hkR2(cIdx,   s) = r2;
        end
        hkConverge(cIdx) = 0;
        fprintf('  HK Case %d  OK\n', caseNum);
    catch ME
        fprintf('  HK Case %d  ERROR: %s\n', caseNum, ME.message);
    end
end

%% Ho-Kalman Aggregate Metrics
hkMeanRMSE  = mean(hkRMSE,  'all', 'omitnan');
hkMeanMAE   = mean(hkMAE,   'all', 'omitnan');
hkMeanR2    = mean(hkR2,    'all', 'omitnan');
hkStability = mean(std(hkRMSE, 0, 1, 'omitnan'), 'omitnan');
hkConvFrac  = mean(hkConverge, 'omitnan');   % algebraic reference (0)

%% Metaheuristic Best-Case Identification
fprintf('\n--- Step 2: Finding best case per algorithm ---\n\n');

results = struct();
results.algoNames    = algorithmNames;
results.bestCaseIdx  = NaN(nAlgos, 1);
results.bestCaseNum  = NaN(nAlgos, 1);
results.RMSE         = NaN(nAlgos, 1);   % 5-state mean at best case
results.Accuracy     = NaN(nAlgos, 1);   % 5-state MAE at best case
results.Stability    = NaN(nAlgos, 1);   % standard deviation of RMSE across cases
results.ConvSpeed    = NaN(nAlgos, 1);   % fractional iteration convergence
results.R2           = NaN(nAlgos, 1);   % 5-state mean R² at best case

algoRMSEAll  = NaN(nAlgos, nCases, 5);
algoR2All    = NaN(nAlgos, nCases, 5);

for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    fprintf('[%d/%d] %s\n', aIdx, nAlgos, algoName);

    % Compute simulation metrics across all available test cases
    for cIdx = 1:nCases
        caseNum = caseNums(cIdx);
        matFile = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, caseNum));
        if ~exist(matFile, 'file'), continue; end
        try
            [X_real, X_sim, ~] = runODE(matFile);
            for s = 1:5
                [rmse, ~, r2] = metrics3(X_real(s,:)', X_sim(:,s));
                algoRMSEAll(aIdx, cIdx, s) = rmse;
                algoR2All(aIdx,  cIdx, s) = r2;
            end
        catch ME
            fprintf('  Case %d  ERROR: %s\n', caseNum, ME.message);
        end
    end

    % Select optimal configuration by highest mean R² across the 5 states
    meanR2PerCase = mean(algoR2All(aIdx, :, :), 3, 'omitnan');
    [~, bestCI]   = max(meanR2PerCase);
    bestCaseNum   = caseNums(bestCI);
    results.bestCaseIdx(aIdx) = bestCI;
    results.bestCaseNum(aIdx) = bestCaseNum;
    fprintf('  → Best case: %s  (mean R² = %.4f)\n', ...
        caseLabels{bestCI}, meanR2PerCase(bestCI));

    % Re-evaluate optimal case trajectory metrics
    bestMatFile = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, bestCaseNum));
    if exist(bestMatFile, 'file')
        try
            [X_real_b, X_sim_b, ~] = runODE(bestMatFile);
            rmseVec = NaN(1,5); maeVec = NaN(1,5); r2Vec = NaN(1,5);
            for s = 1:5
                [rmseVec(s), maeVec(s), r2Vec(s)] = metrics3(X_real_b(s,:)', X_sim_b(:,s));
            end
            results.RMSE(aIdx)     = mean(rmseVec, 'omitnan');
            results.Accuracy(aIdx) = mean(maeVec,  'omitnan');
            results.R2(aIdx)       = mean(r2Vec,   'omitnan');
        catch ME
            fprintf('  Best-case ODE ERROR: %s\n', ME.message);
        end
    end

    % Sensitivity/stability across parameter window variations
    meanRMSEperCase = mean(algoRMSEAll(aIdx, :, :), 3, 'omitnan');
    results.Stability(aIdx) = std(meanRMSEperCase, 'omitnan');

    % Iteration profile convergence index
    results.ConvSpeed(aIdx) = computeConvergence(algoName, bestCaseNum, caseNums);

    fprintf('  RMSE=%.4f  MAE=%.4f  R²=%.4f  Stab=%.4f  Conv=%.4f\n', ...
        results.RMSE(aIdx), results.Accuracy(aIdx), results.R2(aIdx), ...
        results.Stability(aIdx), results.ConvSpeed(aIdx));
end

%% Console Summary Table
fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  BEST CASE COMPARISON: Each Algorithm vs Ho-Kalman\n');
fprintf('  Metrics averaged over 5 longitudinal states\n');
fprintf('================================================================\n\n');

colW = 12;
hdr = sprintf('  %-10s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n', ...
    'Algorithm', 'BestCase', 'Accuracy', 'RMSE', ...
    'Stability', 'ConvSpeed', 'R2', 'Verdict');
fprintf('%s', hdr);
fprintf('  %s\n', repmat('-', 1, 78));

% Ho-Kalman reference baseline
fprintf('  %-10s  %-8s  %-8.4f  %-8.4f  %-8.4f  %-8s  %-8.4f  [Baseline]\n', ...
    'Ho-Kalman', 'N/A', hkMeanMAE, hkMeanRMSE, hkStability, '0 (ref)', hkMeanR2);
fprintf('  %s\n', repmat('-', 1, 78));

for aIdx = 1:nAlgos
    algoName  = algorithmNames{aIdx};
    cLabel    = caseLabels{results.bestCaseIdx(aIdx)};
    acc       = results.Accuracy(aIdx);
    rmse      = results.RMSE(aIdx);
    stab      = results.Stability(aIdx);
    conv      = results.ConvSpeed(aIdx);
    r2        = results.R2(aIdx);

    % Performance categorization based on delta R²
    dR2 = r2 - hkMeanR2;
    if isnan(dR2)
        verdict = 'N/A';
    elseif dR2 > 0.05
        verdict = 'MUCH BETTER';
    elseif dR2 > 0.01
        verdict = 'BETTER';
    elseif dR2 > -0.01
        verdict = 'COMPARABLE';
    elseif dR2 > -0.05
        verdict = 'WORSE';
    else
        verdict = 'MUCH WORSE';
    end

    fprintf('  %-10s  %-8s  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %s\n', ...
        algoName, cLabel, acc, rmse, stab, conv, r2, verdict);
end

fprintf('\n  Notes:\n');
fprintf('  • Accuracy     = Mean Absolute Error (MAE)\n');
fprintf('  • RMSE         = Root Mean Squared Error\n');
fprintf('  • Stability    = Standard deviation of mean RMSE across cases\n');
fprintf('  • ConvSpeed    = Normalized iteration fraction to reach within 1%% of optimum\n');
fprintf('  • R²           = Coefficient of Determination\n');
fprintf('================================================================\n\n');

%% Export Results to Excel
if ~exist('output', 'dir'), mkdir('output'); end
excelOut = fullfile('output', outputExcel);
if exist(excelOut, 'file'), delete(excelOut); end
fprintf('Writing Excel to: %s\n', excelOut);

%% Sheet 1: BestCase_Summary
hdr1 = {'Algorithm','Best_Case', ...
        'Accuracy_MAE','RMSE','Stability_StdRMSE','Conv_Speed_Frac','R2', ...
        'Delta_Accuracy','Delta_RMSE','Delta_Stability','Delta_Conv','Delta_R2', ...
        'Verdict'};
rows1 = {};

% Baseline reference entry
rows1(end+1,:) = {'Ho-Kalman','(all)', ...
    hkMeanMAE, hkMeanRMSE, hkStability, hkConvFrac, hkMeanR2, ...
    0, 0, 0, 0, 0, 'Baseline'};

for aIdx = 1:nAlgos
    algoName = algorithmNames{aIdx};
    cLabel   = caseLabels{results.bestCaseIdx(aIdx)};
    acc  = results.Accuracy(aIdx);
    rmse = results.RMSE(aIdx);
    stab = results.Stability(aIdx);
    conv = results.ConvSpeed(aIdx);
    r2   = results.R2(aIdx);

    dAcc  = acc  - hkMeanMAE;
    dRMSE = rmse - hkMeanRMSE;
    dStab = stab - hkStability;
    dConv = conv - hkConvFrac;
    dR2   = r2   - hkMeanR2;

    dR2v = r2 - hkMeanR2;
    if isnan(dR2v),      verdict = 'N/A';
    elseif dR2v > 0.05,  verdict = 'MUCH BETTER than HK';
    elseif dR2v > 0.01,  verdict = 'BETTER than HK';
    elseif dR2v > -0.01, verdict = 'COMPARABLE to HK';
    elseif dR2v > -0.05, verdict = 'WORSE than HK';
    else,                verdict = 'MUCH WORSE than HK';
    end

    rows1(end+1,:) = {algoName, cLabel, ...
        safeVal(acc), safeVal(rmse), safeVal(stab), safeVal(conv), safeVal(r2), ...
        safeVal(dAcc), safeVal(dRMSE), safeVal(dStab), safeVal(dConv), safeVal(dR2), ...
        verdict}; %#ok<AGROW>
end

T1 = cell2table(rows1, 'VariableNames', hdr1);
writetable(T1, excelOut, 'Sheet', 'BestCase_Summary');
fprintf('  Sheet written: BestCase_Summary\n');

%% Sheet 2: Per-State Detail at Best Case
stateNames = {'u(m/s)','w(m/s)','q(rad/s)','theta(rad)','h(m)'};
hdr2 = [{'Algorithm','Best_Case','State'}, ...
        {'HK_RMSE','HK_MAE','HK_R2'}, ...
        {'Algo_RMSE','Algo_MAE','Algo_R2'}, ...
        {'Delta_RMSE','Delta_MAE','Delta_R2','Winner'}];
rows2 = {};

for aIdx = 1:nAlgos
    algoName    = algorithmNames{aIdx};
    bestCaseNum = results.bestCaseNum(aIdx);
    bestCI      = results.bestCaseIdx(aIdx);
    cLabel      = caseLabels{bestCI};
    bestMatFile = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, bestCaseNum));
    hkFile      = fullfile('output', sprintf('HOKALMAN_Case%d_24.mat', bestCaseNum));

    if ~exist(bestMatFile,'file') || ~exist(hkFile,'file'), continue; end
    try
        [Xr_a, Xs_a, ~] = runODE(bestMatFile);
        [Xr_h, Xs_h, ~] = runODE(hkFile);
        for s = 1:5
            [ar, am, aR2] = metrics3(Xr_a(s,:)', Xs_a(:,s));
            [hr, hm, hR2] = metrics3(Xr_h(s,:)', Xs_h(:,s));
            dR2w = aR2 - hR2;
            if isnan(dR2w),     win = 'N/A';
            elseif dR2w > 0.001, win = algoName;
            elseif dR2w < -0.001,win = 'HOKALMAN';
            else,                win = 'Draw';
            end
            rows2(end+1,:) = {algoName, cLabel, stateNames{s}, ...
                hr, hm, hR2, ar, am, aR2, ...
                ar-hr, am-hm, dR2w, win}; %#ok<AGROW>
        end
    catch ME
        fprintf('  Per-state detail ERROR (%s): %s\n', algoName, ME.message);
    end
end
T2 = cell2table(rows2, 'VariableNames', hdr2);
writetable(T2, excelOut, 'Sheet', 'PerState_BestCase');
fprintf('  Sheet written: PerState_BestCase\n');

%% Sheet 3: Convergence History
hdr3 = [{'Algorithm','Best_Case','MaxIter','Conv_Iter','Conv_Frac'}, ...
         arrayfun(@(i) sprintf('Iter_%d', i), 1:20, 'UniformOutput', false)];
rows3 = {};
for aIdx = 1:nAlgos
    algoName    = algorithmNames{aIdx};
    bestCI      = results.bestCaseIdx(aIdx);
    bestCaseNum = results.bestCaseNum(aIdx);
    cLabel      = caseLabels{bestCI};

    matFile = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, bestCaseNum));
    if ~exist(matFile,'file'), continue; end
    try
        S = load(matFile, 'globalFitnessMatrix', 'maxIterations', 'iter');
        if ~isfield(S,'globalFitnessMatrix'), continue; end

        if isfield(S,'iter'), nIt = S.iter;
        elseif isfield(S,'maxIterations'), nIt = S.maxIterations;
        else, nIt = size(S.globalFitnessMatrix,1); end

        meanFit = zeros(nIt,1);
        for it = 1:nIt
            meanFit(it) = mean(S.globalFitnessMatrix(it,:), 'omitnan');
        end
        [convIter, convFrac] = convergenceIndex(meanFit, nIt);

        padded = NaN(1,20);
        padded(1:min(nIt,20)) = meanFit(1:min(nIt,20))';

        rows3(end+1,:) = [{algoName, cLabel, nIt, convIter, convFrac}, num2cell(padded)]; %#ok<AGROW>
    catch ME
        fprintf('  Convergence detail ERROR (%s): %s\n', algoName, ME.message);
    end
end
if ~isempty(rows3)
    T3 = cell2table(rows3, 'VariableNames', hdr3);
    writetable(T3, excelOut, 'Sheet', 'Convergence_Detail');
    fprintf('  Sheet written: Convergence_Detail\n');
end

%% Sheet 4: Composite Ranking
hdr4 = {'Rank','Algorithm','Best_Case', ...
        'R2','RMSE','Accuracy_MAE','Stability','ConvSpeed', ...
        'R2_Rank','RMSE_Rank','Acc_Rank','Stab_Rank','Conv_Rank', ...
        'Composite_Rank','Verdict'};

allNames = [algorithmNames, {'Ho-Kalman'}];
allR2    = [results.R2;    hkMeanR2];
allRMSE  = [results.RMSE;  hkMeanRMSE];
allAcc   = [results.Accuracy; hkMeanMAE];
allStab  = [results.Stability; hkStability];
allConv  = [results.ConvSpeed; hkConvFrac];
bestCase_labels = [caseLabels(results.bestCaseIdx)]; 
bestCase_labels{end+1} = '(all)';

nTotal = numel(allNames);

[~, oR2]   = sort(allR2,   'descend', 'MissingPlacement','last');
[~, oRMSE] = sort(allRMSE, 'ascend',  'MissingPlacement','last');
[~, oAcc]  = sort(allAcc,  'ascend',  'MissingPlacement','last');
[~, oStab] = sort(allStab, 'ascend',  'MissingPlacement','last');
[~, oConv] = sort(allConv, 'ascend',  'MissingPlacement','last');

rankR2   = NaN(nTotal,1); rankR2(oR2)     = 1:nTotal;
rankRMSE = NaN(nTotal,1); rankRMSE(oRMSE) = 1:nTotal;
rankAcc  = NaN(nTotal,1); rankAcc(oAcc)   = 1:nTotal;
rankStab = NaN(nTotal,1); rankStab(oStab) = 1:nTotal;
rankConv = NaN(nTotal,1); rankConv(oConv) = 1:nTotal;
compositeRank = mean([rankR2, rankRMSE, rankAcc, rankStab, rankConv], 2, 'omitnan');
[~, finalOrd] = sort(compositeRank);

rows4 = {};
for r = 1:nTotal
    i = finalOrd(r);
    dR2v = allR2(i) - hkMeanR2;
    if i == nTotal
        verdict = 'Baseline';
    elseif isnan(dR2v),      verdict = 'N/A';
    elseif dR2v > 0.05,  verdict = 'MUCH BETTER than HK';
    elseif dR2v > 0.01,  verdict = 'BETTER than HK';
    elseif dR2v > -0.01, verdict = 'COMPARABLE to HK';
    elseif dR2v > -0.05, verdict = 'WORSE than HK';
    else,                verdict = 'MUCH WORSE than HK';
    end
    rows4(end+1,:) = {r, allNames{i}, bestCase_labels{i}, ...
        safeVal(allR2(i)), safeVal(allRMSE(i)), safeVal(allAcc(i)), ...
        safeVal(allStab(i)), safeVal(allConv(i)), ...
        rankR2(i), rankRMSE(i), rankAcc(i), rankStab(i), rankConv(i), ...
        compositeRank(i), verdict}; %#ok<AGROW>
end
T4 = cell2table(rows4, 'VariableNames', hdr4);
writetable(T4, excelOut, 'Sheet', 'Overall_Ranking');
fprintf('  Sheet written: Overall_Ranking\n');

fprintf('\n================================================================\n');
fprintf('  DONE.  Results saved to:\n  %s\n', fullfile(pwd,'output',outputExcel));
fprintf('================================================================\n\n');

end


%% Helper Functions

function [X_real, X_sim, time] = runODE(matFile)
% Simulates longitudinal flight dynamics via ode15s across flight test window
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

    % Extract longitudinal states: [u, w, q, theta, h]
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
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

    [t_sim, X_raw] = ode15s(odefun, time, X0, opts);
    X_sim = interp1(t_sim, X_raw, time, 'linear', 'extrap');
end

function dx = LongODE_tv_local(t, x, A, B, time, u, elevator, flapPos, flapNeg, flapDiff, thrust)
% Longitudinal state equations with time-varying u(t) for kinematic altitude rate
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

function [rmse, mae, r2] = metrics3(y, yhat)
% Computes goodness-of-fit metrics (RMSE, MAE, R²) for state vectors
    y    = y(:); yhat = yhat(:);
    valid = isfinite(y) & isfinite(yhat);
    y = y(valid); yhat = yhat(valid);
    if isempty(y)
        rmse = NaN; mae = NaN; r2 = NaN;
        return;
    end
    err    = y - yhat;
    rmse   = sqrt(mean(err.^2));
    mae    = mean(abs(err));
    SS_res = sum(err.^2);
    SS_tot = sum((y - mean(y)).^2);
    if SS_tot == 0
        r2 = NaN;
    else
        r2 = 1 - SS_res / SS_tot;
    end
end

function convFrac = computeConvergence(algoName, bestCaseNum, caseNums)
% Extracts convergence fraction from optimization history
    convFrac = NaN;
    matFile  = fullfile('output', sprintf('%s_Case%d_24.mat', algoName, bestCaseNum));
    if ~exist(matFile, 'file'), return; end
    try
        S = load(matFile, 'globalFitnessMatrix', 'maxIterations', 'iter');
        if ~isfield(S,'globalFitnessMatrix'), return; end

        if isfield(S,'iter'), nIt = S.iter;
        elseif isfield(S,'maxIterations'), nIt = S.maxIterations;
        else, nIt = size(S.globalFitnessMatrix, 1); end

        meanFit = zeros(nIt, 1);
        for it = 1:nIt
            meanFit(it) = mean(S.globalFitnessMatrix(it,:), 'omitnan');
        end
        [~, convFrac] = convergenceIndex(meanFit, nIt);
    catch
        convFrac = NaN;
    end
end

function [convIter, convFrac] = convergenceIndex(meanFit, nIt)
% Identifies iteration at which mean fitness reaches within 1% of final value
    finalVal = meanFit(end);
    initVal  = meanFit(1);
    threshold = finalVal + 0.01 * abs(initVal - finalVal);
    convIter  = nIt;
    for it = 1:nIt
        if meanFit(it) <= threshold
            convIter = it;
            break;
        end
    end
    convFrac = convIter / nIt;
end

function v = safeVal(x)
% Sanitizes non-finite values for Excel cell writing
    if isnan(x) || isinf(x)
        v = 'N/A';
    else
        v = x;
    end
end
