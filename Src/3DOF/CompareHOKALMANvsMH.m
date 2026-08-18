function CompareHOKALMANvsMH(A_HK, B_HK, A_MH, B_MH, varargin)
% CompareHOKALMANvsMH Compares time-domain state validation of Ho-Kalman vs
% metaheuristic identified state-space matrices against recorded flight data.
%
% Inputs:
%   A_HK, B_HK  - Longitudinal system matrices (5x5) from Ho-Kalman identification
%   A_MH, B_MH  - Longitudinal system matrices (5x5) from metaheuristic optimization
%
% Optional Name-Value Pairs:
%   'CaseNumber'       - Case identifier index (default: 11)
%   'AlgorithmName'    - Name of metaheuristic optimizer (default: 'MH')
%   'OutputDir'        - Destination folder for exported plots (default: 'output')
%   'ValidationMargin' - Boundary trim points to exclude transient edges (default: 1000)

%% Parse Arguments
p = inputParser;
addRequired(p, 'A_HK', @(x) isequal(size(x), [5 5]));
addRequired(p, 'B_HK', @(x) isequal(size(x), [5 5]));
addRequired(p, 'A_MH', @(x) isequal(size(x), [5 5]));
addRequired(p, 'B_MH', @(x) isequal(size(x), [5 5]));
addParameter(p, 'CaseNumber', 11, @isscalar);
addParameter(p, 'AlgorithmName', 'MH', @ischar);
addParameter(p, 'OutputDir', 'output', @ischar);
addParameter(p, 'ValidationMargin', 1000, @isscalar);
parse(p, A_HK, B_HK, A_MH, B_MH, varargin{:});

caseNumber = p.Results.CaseNumber;
algoName = p.Results.AlgorithmName;
outDir = p.Results.OutputDir;
margin = p.Results.ValidationMargin;

%% Load Flight Data
accelFile = "cruise_acceleration.mat";
statesFile = "cruise_outputStates.mat";
accutFile = "cruise_acctuators.mat";
accel = load(accelFile);
states = load(statesFile);
accut = load(accutFile);

popNum = 10100; offset = 2000; limit = popNum + offset;
pdot_qdot_rdot = states.output_states.pdot_qdot_rdot.Data(offset:limit, :);
pqr = states.output_states.pqr.Data(offset:limit, :);
phi_theta_psi = squeeze(states.output_states.phi_theta_psi.Data);
phi_theta_psi = phi_theta_psi(:, offset:limit);
Vb = states.output_states.Vb.Data(offset:limit, :);
Ve = states.output_states.Ve.Data(offset:limit, :);
Xe = states.output_states.Xe.Data(offset:limit, :);
Accels = accel.acceleration_bb.Data(offset:limit, :);
elevator = accut.accutators.Data(offset:limit, 2);
time = states.output_states.pdot_qdot_rdot.Time(offset:limit, :) - ...
       states.output_states.pdot_qdot_rdot.Time(offset, :);

%% Validation Margin Trimming
N = size(Vb, 1);
if margin > 0
    idxS = 1 + margin; idxE = N - margin;
    Vb = Vb(idxS:idxE, :); pqr = pqr(idxS:idxE, :);
    phi_theta_psi = phi_theta_psi(:, idxS:idxE);
    Xe = Xe(idxS:idxE, :); Accels = Accels(idxS:idxE, :);
    elevator = elevator(idxS:idxE, :); Ve = Ve(idxS:idxE, :);
    pdot_qdot_rdot = pdot_qdot_rdot(idxS:idxE, :);
    time = time(idxS:idxE); time = time - time(1);
end

% Assemble state vector [u, w, q, theta, h]
N = size(Vb, 1);
u_data = Vb(:, 1); w_data = Vb(:, 3);
q_data = pqr(:, 2); theta_data = phi_theta_psi(2, :)';
h_data = -Xe(:, 3);
elevator_data = elevator(:, 1);
X_real = [u_data, w_data, q_data, theta_data, h_data]';
X0 = X_real(:, 1);

thrust = 0;
state_names = {'$u$ (m/s)', '$w$ (m/s)', '$q$ (rad/s)', '$\theta$ (rad)', '$h$ (m)'};
dt = mean(diff(time));

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

caseDisplay = sprintf('(%d,%d)', floor(caseNumber/10), mod(caseNumber,10));

%% Numerical ODE Simulation (ode15s)
odefun_HK = @(t, x) LongODE_tv(t, x, A_HK, B_HK, time, u_data, elevator_data, thrust);
odefun_MH = @(t, x) LongODE_tv(t, x, A_MH, B_MH, time, u_data, elevator_data, thrust);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

[t_HK, X_HK_raw] = ode15s(odefun_HK, time, X0, options);
X_HK = interp1(t_HK, X_HK_raw, time);
[t_MH, X_MH_raw] = ode15s(odefun_MH, time, X0, options);
X_MH = interp1(t_MH, X_MH_raw, time);

%% Error Metrics
fprintf('\n=========================================================\n');
fprintf('  Comparison: Ho-Kalman vs %s - Case %s\n', algoName, caseDisplay);
fprintf('=========================================================\n');
fprintf('\n--- ODE State Simulation RMSE ---\n');
fprintf('  %-12s  %14s  %14s\n', 'State', 'Ho-Kalman', algoName);
fprintf('  %s\n', repmat('-', 1, 44));
for i = 1:5
    rmse_HK = sqrt(mean((X_real(i,:)' - X_HK(:,i)).^2));
    rmse_MH = sqrt(mean((X_real(i,:)' - X_MH(:,i)).^2));
    fprintf('  %-12s  %14.4f  %14.4f\n', state_names{i}(2:end-4), rmse_HK, rmse_MH);
end

%% Plot 1: State Trajectory Comparison
fig1 = figure('Units', 'inches', 'Position', [1 1 7.16 7.5], 'Color', 'w');
for i = 1:5
    subplot(3, 2, i);
    plot(time, X_real(i, :), 'k-', 'LineWidth', 0.8); hold on;
    plot(time, X_HK(:, i), '--b', 'LineWidth', 0.8);
    plot(time, X_MH(:, i), '--r', 'LineWidth', 0.8);
    xlabel('Time (s)'); ylabel(state_names{i}, 'Interpreter', 'latex');
    legend({'Flight', 'Ho-Kalman', algoName}, 'Location', 'best', 'FontSize', 7);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
end
sgtitle(sprintf('State Validation: Ho-Kalman vs %s - Case %s', algoName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig1, fullfile(outDir, sprintf('Compare_States_HKvs%s_Case%d_Longitudinal.png', algoName, caseNumber)), 'Resolution', 600);

%% Plot 2: Residual Error Comparison
fig2 = figure('Units', 'inches', 'Position', [1 1 7.16 7.5], 'Color', 'w');
for i = 1:5
    subplot(3, 2, i);
    err_HK = X_real(i, :)' - X_HK(:, i);
    err_MH = X_real(i, :)' - X_MH(:, i);
    plot(time, err_HK, 'b', 'LineWidth', 0.6); hold on;
    plot(time, err_MH, 'r', 'LineWidth', 0.6);
    xlabel('Time (s)'); ylabel('Error');
    title(state_names{i}, 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 9);
    legend({'Ho-Kalman', algoName}, 'Location', 'best', 'FontSize', 7);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
end
sgtitle(sprintf('Validation Error: Ho-Kalman vs %s - Case %s', algoName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig2, fullfile(outDir, sprintf('Compare_Error_HKvs%s_Case%d_Longitudinal.png', algoName, caseNumber)), 'Resolution', 600);

fprintf('\nFigures saved to %s/\n', outDir);
fprintf('  Compare_States_HKvs%s_Case%d_Longitudinal.png\n', algoName, caseNumber);
fprintf('  Compare_Error_HKvs%s_Case%d_Longitudinal.png\n', algoName, caseNumber);
fprintf('\nDone.\n');

end

%% Local Helper: Longitudinal Equations with Time-Varying Speed
function dx = LongODE_tv(t, x, A, B, time, u, elevator, thrust)
    de = interp1(time, elevator, t, 'linear', 'extrap');
    ut = interp1(time, u, t, 'linear', 'extrap');
    U = [de; thrust; 0; 0; 0];
    A_tv = A; A_tv(5, 4) = ut;
    dx = A_tv * x + B * U;
end
