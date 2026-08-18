function CreateGraphicsLongitudinalDynamics24(outputFileName, caseNumber, algorithmName)
% CreateGraphicsLongitudinalDynamics24 Generates publication-ready validation
% plots and metrics for the 24-parameter longitudinal aircraft model, including
% ODE15s state simulation, one-step derivative checks, modal eigenvalues,
% step/impulse responses, and Excel report generation.
%
% Inputs:
%   outputFileName - Saved MAT workspace file from identification run
%   caseNumber     - Test matrix case identifier
%   algorithmName  - Optimization algorithm name (e.g., 'ALO')

%% Workspace and Environment Setup
load(outputFileName);

try
    opengl software;
catch
end

caseDisplay = sprintf('(%d,%d)', floor(caseNumber/10), mod(caseNumber,10));

%% Figure Formatting Defaults
set(0, 'DefaultFigureVisible', 'off');
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
set(0, 'DefaultLegendFontSize', 8);
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');

g = valueOfGravitationConstant;

try
    margin = inputData.ValidationMargin;
catch
    margin = 0;
end

N = size(Vb, 1);
idxStart = 1 + margin;
idxEnd   = N - margin;
Nvalid   = idxEnd - idxStart + 1;

Vtrim = Vb(idxStart:idxEnd, 1);
uo_mean = mean(Vtrim);

thrust = 0;

[A_id, B_id] = formatParameters24(globalBest, uo_mean);

%% Flight Data Extraction and Trimming
if isfield(staticData, 'time_all')
    time = staticData.time_all;
    dataIsTrimmed = true;
else
    time = (0:N-1)';
    dataIsTrimmed = false;
end

dStart = 1 + margin;
dEnd   = N - margin;

u_data    = Vb(   dStart:dEnd, 1);
w_data    = Vb(   dStart:dEnd, 3);
q_data    = pqr(  dStart:dEnd, 2);
theta_data = phi_theta_psi(2, dStart:dEnd)';
h_data    = -Xe(  dStart:dEnd, 3);
elevator_data = elevator(dStart:dEnd, 1);
flapPos_data  = flapPos( dStart:dEnd, 1);
flapNeg_data  = flapNeg( dStart:dEnd, 1);
flapDiff_data = flapDiff(dStart:dEnd, 1);

if dataIsTrimmed
    time = time;
else
    time = time(dStart:dEnd);
end
dt = mean(diff(time));

% Longitudinal state vector [u, w, q, theta, h]
X_real = [u_data, w_data, q_data, theta_data, h_data]';
X0 = X_real(:, 1);

state_names = {'$u$ (m/s)', '$w$ (m/s)', '$q$ (rad/s)', '$\theta$ (rad)', '$h$ (m)'};

%% ODE15s State Simulation
odefun = @(t, x) LongODE_tv(t, x, A_id, B_id, time, u_data, elevator_data, flapPos_data, flapNeg_data, flapDiff_data, thrust);

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t_sim, X_sim_raw] = ode15s(odefun, time, X0, options);
X_sim = interp1(t_sim, X_sim_raw, time);

%% One-Step Derivative Prediction
Nobs = size(X_real, 2);
X_dot_model = A_id * X_real + B_id * [elevator_data'; thrust * ones(1, Nobs); flapPos_data'; flapNeg_data'; flapDiff_data'];

if dataIsTrimmed
    x_dot_true = staticData.sim_error_check_all;
else
    x_dot_true = staticData.sim_error_check_all(:, dStart:dEnd);
end

%% Metric Calculations and Console Report
fprintf('\n====================================\n');
fprintf('%s - Case %s - VALIDATION REPORT\n', algorithmName, caseDisplay);
fprintf('====================================\n');

fprintf('\n--- ODE State Simulation (A matrix with uo = mean u) ---\n');
MetricTableODE = zeros(5, 4);
for i = 1:5
    y = X_real(i, :)';
    yhat = X_sim(:, i);
    RMSE = sqrt(mean((y - yhat).^2));
    MAE = mean(abs(y - yhat));
    R2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    FIT = 100 * (1 - norm(y - yhat) / norm(y - mean(y)));
    fprintf('  %s  RMSE=%.4f  MAE=%.4f  R2=%.4f  FIT=%.2f%%\n', state_names{i}, RMSE, MAE, R2, FIT);
    MetricTableODE(i, :) = [RMSE, MAE, R2, FIT];
end

fprintf('\n--- One-Step Derivative Prediction (time-varying u, matches training) ---\n');
MetricTable1step = zeros(5, 4);
for i = 1:5
    y = x_dot_true(i, :)';
    yhat = X_dot_model(i, :)';
    RMSE = sqrt(mean((y - yhat).^2));
    MAE = mean(abs(y - yhat));
    R2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    FIT = 100 * (1 - norm(y - yhat) / norm(y - mean(y)));
    fprintf('  dot_%s  RMSE=%.4f  MAE=%.4f  R2=%.4f  FIT=%.2f%%\n', ...
        state_names{i}(2), RMSE, MAE, R2, FIT);
    MetricTable1step(i, :) = [RMSE, MAE, R2, FIT];
end

%% Plot 1: RMSE Convergence by Iteration Block
fig = figure('Units', 'inches', 'Position', [1 1 7.16 5.5], 'Color', 'w');
iVec = 1:populationSize;
titles = {'Iterations 1-5', 'Iterations 6-10', 'Iterations 11-15', 'Iterations 16-20'};
itBlocks = {1:5, 6:10, 11:15, 16:20};

for subIdx = 1:4
    subplot(2, 2, subIdx);
    j = itBlocks{subIdx};
    j = j(j <= iter);
    if ~isempty(j)
        plot(iVec, globalFitnessMatrix(j, :)', 'LineWidth', 0.6);
        legend(cellstr(num2str(j', 'it-%d')), 'Location', 'northeast', 'FontSize', 7);
    end
    xlim([1 populationSize]);
    ylim([0 1500]);
    xlabel('Population Index');
    ylabel('RMSE');
    title(titles{subIdx}, 'FontWeight', 'normal', 'FontSize', 9);
end

han = axes(fig, 'visible', 'off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, 'RMSE');
xlabel(han, 'Population Index');
sgtitle(sprintf('%s - Case %s - RMSE Convergence', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig, sprintf('output/Error%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig, sprintf('output/Error%s_Case%d_Longitudinal.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Plot 2: Time-Domain State Validation
fig2 = figure('Units', 'inches', 'Position', [1 1 7.16 6.5], 'Color', 'w');
for i = 1:5
    subplot(3, 2, i);
    plot(time, X_real(i, :), 'k-', 'LineWidth', 1.2); hold on;
    plot(time, X_sim(:, i), '--r', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel(state_names{i}, 'Interpreter', 'latex');
    legend({'Flight', 'Simulated'}, 'Location', 'best', 'FontSize', 7);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
end
sgtitle(sprintf('%s - Case %s - State Validation', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig2, sprintf('output/StateValidation_%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);

%% Plot 3: Residual Errors
fig3 = figure('Units', 'inches', 'Position', [1 1 7.16 6.5], 'Color', 'w');
for i = 1:5
    subplot(3, 2, i);
    err = X_real(i, :)' - X_sim(:, i);
    plot(time, err, 'k', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Error');
    title(state_names{i}, 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 9);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
end
sgtitle(sprintf('%s - Case %s - Validation Error', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig3, sprintf('output/ErrorResiduals_%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);

%% Plot 4: System Modal Eigenvalues
lambda = eig(A_id);
fprintf('\n========== Eigenvalues ==========\n');
disp(lambda);

fig4 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 3.5 2.8]);
plot(real(lambda), imag(lambda), 'x', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Identified Aircraft Modes', 'FontWeight', 'normal', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
exportgraphics(fig4, sprintf('output/Eigenvalues_%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);

%% Plot 5: Step Response
sys_id = ss(A_id, B_id, eye(5), zeros(5, 5));
sys_uwq = sys_id(1:3, :);
fig5 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7.16 5.5]);
step(sys_uwq);
set(findall(gcf, 'Type', 'Axes'), 'FontName', 'Times New Roman', 'FontSize', 9);
title(sprintf('%s - Case %s - Step Response (u, w, q)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig5, sprintf('output/StepResponse_%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);

%% Plot 6: Impulse Response
fig6 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7.16 5.5]);
impulse(sys_uwq);
set(findall(gcf, 'Type', 'Axes'), 'FontName', 'Times New Roman', 'FontSize', 9);
title(sprintf('%s - Case %s - Impulse Response (u, w, q)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig6, sprintf('output/ImpulseResponse_%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);

%% Plot 7: Mean Learning Curve
fig7 = figure('Units', 'inches', 'Position', [1 1 3.5 2.8], 'Color', 'w');
avgRMSE = zeros(iter, 1);
for k = 1:iter
    avgRMSE(k) = mean(globalFitnessMatrix(k, :));
end
plot(1:iter, avgRMSE, 'k-', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('Mean RMSE');
title('Mean RMSE vs. Iteration', 'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig7, sprintf('output/MeanErrorGraph%s_Case%d_Longitudinal.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig7, sprintf('output/MeanErrorGraph%s_Case%d_Longitudinal.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Export Validation Metrics to Excel
excelFile = sprintf('output/ValidationReport_%s_Case%d.xlsx', algorithmName, caseNumber);
ResultsODE = table(MetricTableODE(:, 1), MetricTableODE(:, 2), MetricTableODE(:, 3), MetricTableODE(:, 4), ...
    'VariableNames', {'RMSE', 'MAE', 'R2', 'FIT_Percent'}, ...
    'RowNames', {'u', 'w', 'q', 'theta', 'h'});
Results1step = table(MetricTable1step(:, 1), MetricTable1step(:, 2), MetricTable1step(:, 3), MetricTable1step(:, 4), ...
    'VariableNames', {'RMSE', 'MAE', 'R2', 'FIT_Percent'}, ...
    'RowNames', {'u_dot', 'w_dot', 'q_dot', 'theta_dot', 'h_dot'});
writetable(ResultsODE, excelFile, 'Sheet', 'ODE_State_Simulation', 'WriteRowNames', true);
writetable(Results1step, excelFile, 'Sheet', 'OneStep_Derivative', 'WriteRowNames', true);

%% Export Optimal Parameters and Friedman Test
bestSheet = sprintf('output/Data_Case%d_Longitudinal.xlsx', caseNumber);
try
    xlswrite(bestSheet, globalBest, 'Sheet1');
    gf = globalFitnessMatrix';
    [p, tbl, stats] = friedman(gf, iter, 'off');
    tblArray = cell2table(tbl(2:end, :), 'VariableNames', tbl(1, :));
    writetable(tblArray, bestSheet, 'Sheet', 'Sheet2');
catch ME
    warning('Failed Friedman test: %s', ME.message);
    try; xlswrite(bestSheet, globalBest, 'Sheet1'); catch; end
end

fprintf('\n========== VALIDATION COMPLETE ==========\n');
fprintf('All figures and report saved to output/\n');

end

%% Helper Functions: Longitudinal ODE Equations

function dx = LongODE(t, x, A, B, time, elevator, flapPos, flapNeg, flapDiff, thrust)
% Linear longitudinal state equation with constant system matrices
    de = interp1(time, elevator, t, 'linear', 'extrap');
    fp = interp1(time, flapPos, t, 'linear', 'extrap');
    fn = interp1(time, flapNeg, t, 'linear', 'extrap');
    fd = interp1(time, flapDiff, t, 'linear', 'extrap');
    U = [de; thrust; fp; fn; fd];
    dx = A * x + B * U;
end

function dx = LongODE_tv(t, x, A, B, time, u, elevator, flapPos, flapNeg, flapDiff, thrust)
% Longitudinal state equation with time-varying u(t) coupling in kinematic altitude rate
    de = interp1(time, elevator, t, 'linear', 'extrap');
    fp = interp1(time, flapPos, t, 'linear', 'extrap');
    fn = interp1(time, flapNeg, t, 'linear', 'extrap');
    fd = interp1(time, flapDiff, t, 'linear', 'extrap');
    ut = interp1(time, u, t, 'linear', 'extrap');
    U = [de; thrust; fp; fn; fd];
    A_tv = A;
    A_tv(5, 4) = ut;
    dx = A_tv * x + B * U;
end
