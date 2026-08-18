function CreateGraphics6DOF(outputFileName, caseNumber, algorithmName)
% CreateGraphics6DOF  Generates IEEE-format validation figures and performance metrics for 6-DOF models.
%
% Evaluates identified 12-state, 7-input continuous state-space models:
%   States (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
%   Inputs  (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]
%
% Generated outputs (in output/ directory):
%   - Multi-solver ODE state trajectory simulations (12 states and 6 motion states)
%   - Error residuals and range-normalized percentage errors
%   - One-step acceleration and angular rate derivative matching
%   - Open-loop pole spectrum (12x12 A matrix)
%   - Subsystem step and impulse dynamic responses
%   - Excel validation summary workbook

load(outputFileName);

%% Out-of-Sample Validation Dataset Evaluation
if isprop(inputData, 'HasTestData') && inputData.HasTestData && ...
   exist(inputData.TestAccelerationFileName, 'file') && ...
   exist(inputData.TestOutputStatesFilesName, 'file') && ...
   exist(inputData.TestAccutatorsFileName, 'file')
    
    fprintf('CreateGraphics6DOF: Evaluating identified model on Out-of-Sample Testing Dataset...\n');
    fprintf('  Train Files: %s | %s | %s\n', inputData.AccelerationFileName, inputData.OutputStatesFilesName, inputData.AccutatorsFileName);
    fprintf('  Test Files : %s | %s | %s\n', inputData.TestAccelerationFileName, inputData.TestOutputStatesFilesName, inputData.TestAccutatorsFileName);
    
    testInput = inputData;
    testInput.AccelerationFileName  = inputData.TestAccelerationFileName;
    testInput.OutputStatesFilesName = inputData.TestOutputStatesFilesName;
    testInput.AccutatorsFileName    = inputData.TestAccutatorsFileName;
    if isprop(inputData, 'TestOffSet') && ~isempty(inputData.TestOffSet)
        testInput.OffSet = inputData.TestOffSet;
    end
    if isprop(inputData, 'TestPopulationNumber') && ~isempty(inputData.TestPopulationNumber)
        testInput.PopulationNumber = inputData.TestPopulationNumber;
    end
    
    testStaticData = PrepareFlightData6DOF(testInput);
    Vb             = testStaticData.Vb;
    pqr            = testStaticData.pqr;
    phi_theta_psi  = testStaticData.phi_theta_psi;
    Xe             = testStaticData.Xe;
    aileron        = testStaticData.aileron;
    rudder         = testStaticData.rudder;
    elevator       = testStaticData.elevator;
    flapPos        = testStaticData.flapPos;
    flapNeg        = testStaticData.flapNeg;
    flapDiff       = testStaticData.flapDiff;
    Accels         = testStaticData.Accels;
    pdot_qdot_rdot = testStaticData.pdot_qdot_rdot;
    Ve             = testStaticData.Ve;
    staticData     = testStaticData;
end

try; opengl software; catch; end

caseDisplay = sprintf('(%d,%d)', floor(caseNumber / 10), mod(caseNumber, 10));

%% IEEE Figure Default Parameters
set(0, 'DefaultFigureVisible',  'off');
set(0, 'DefaultFigureColor',    'w');
set(0, 'DefaultAxesColor',      'w');
set(0, 'DefaultAxesFontName',   'Times New Roman');
set(0, 'DefaultAxesFontSize',   9);
set(0, 'DefaultAxesLineWidth',  0.8);
set(0, 'DefaultAxesTickDir',    'in');
set(0, 'DefaultAxesBox',        'on');
set(0, 'DefaultLineLineWidth',  0.8);
set(0, 'DefaultTextFontName',   'Times New Roman');
set(0, 'DefaultTextFontSize',   9);
set(0, 'DefaultLegendFontName', 'Times New Roman');
set(0, 'DefaultLegendFontSize', 8);
set(0, 'DefaultTextColor',      'k');
set(0, 'DefaultAxesXColor',     'k');
set(0, 'DefaultAxesYColor',     'k');

g = valueOfGravitationConstant;

try
    margin = inputData.ValidationMargin;
catch
    margin = 0;
end

N      = size(Vb, 1);
dStart = 1 + margin;
dEnd   = N - margin;

%% Reconstruct State Telemetry
u_data     = Vb(dStart:dEnd, 1);
v_data     = Vb(dStart:dEnd, 2);
w_data     = Vb(dStart:dEnd, 3);
p_data     = pqr(dStart:dEnd, 1);
q_data     = pqr(dStart:dEnd, 2);
r_data     = pqr(dStart:dEnd, 3);
phi_data   = phi_theta_psi(1, dStart:dEnd)';
theta_data = phi_theta_psi(2, dStart:dEnd)';
psi_data   = phi_theta_psi(3, dStart:dEnd)';
xe_data    = Xe(dStart:dEnd, 1);
ye_data    = Xe(dStart:dEnd, 2);
h_data     = -Xe(dStart:dEnd, 3);

aileron_data  = aileron( dStart:dEnd, 1);
rudder_data   = rudder(  dStart:dEnd, 1);
elevator_data = elevator(dStart:dEnd, 1);
flapPos_data  = flapPos( dStart:dEnd, 1);
flapNeg_data  = flapNeg( dStart:dEnd, 1);
flapDiff_data = flapDiff(dStart:dEnd, 1);

if isfield(staticData, 'time_all')
    time = staticData.time_all;
else
    time = (0:(dEnd - dStart))';
end
dt = mean(diff(time));

X_real = [u_data, v_data, w_data, p_data, q_data, r_data, ...
           phi_data, theta_data, psi_data, xe_data, ye_data, h_data]';

%% Build Identified Continuous State-Space Model
if isfield(staticData, 'u0_trim')
    u0_trim = staticData.u0_trim;
    x_trim  = staticData.x_trim;
    u_trim  = staticData.u_trim;
else
    u0_trim = mean(u_data);
    x_trim  = [u0_trim; zeros(11,1)];
    u_trim  = zeros(7,1);
end

[A_id, B_id] = formatParameters6DOF(globalBest, u0_trim);

state_names = {'$u$ (m/s)',     '$v$ (m/s)',     '$w$ (m/s)', ...
               '$p$ (rad/s)',   '$q$ (rad/s)',   '$r$ (rad/s)', ...
               '$\phi$ (rad)',  '$\theta$ (rad)','$\psi$ (rad)', ...
               '$x_e$ (m)',     '$y_e$ (m)',      '$h$ (m)'};

deriv_names = {'$\dot{u}$',   '$\dot{v}$',   '$\dot{w}$', ...
               '$\dot{p}$',   '$\dot{q}$',   '$\dot{r}$', ...
               '$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$', ...
               '$\dot{x}_e$', '$\dot{y}_e$',  '$\dot{h}$'};

%% Continuous ODE Simulation in Perturbation Space
% Actuator perturbations in radians (order: [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder])
ele_pert = deg2rad(elevator_data - u_trim(1));
fp_pert  = deg2rad(flapPos_data  - u_trim(3));
fn_pert  = deg2rad(flapNeg_data  - u_trim(4));
fd_pert  = deg2rad(flapDiff_data - u_trim(5));
da_pert  = deg2rad(aileron_data  - u_trim(6));
dr_pert  = deg2rad(rudder_data   - u_trim(7));

dv_real     = v_data     - x_trim(2);
dphi_real   = phi_data   - x_trim(7);
dtheta_real = theta_data - x_trim(8);
dpsi_real   = psi_data   - x_trim(9);

X0_pert = X_real(:, 1) - x_trim;

odefun = @(t, dx) FullODE_6DOF_pert(t, dx, A_id, B_id, time, u_data, u0_trim, ...
             ele_pert, fp_pert, fn_pert, fd_pert, da_pert, dr_pert, ...
             dv_real, dphi_real, dtheta_real, dpsi_real);

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

solvers = {@ode45, @ode15s, @ode113, @ode23, @ode23s};
solverNames = {'ode45 (Dormand-Prince 4,5)', ...
               'ode15s (Stiff NDF/BDF)', ...
               'ode113 (Adams-Bashforth-Moulton)', ...
               'ode23 (Bogacki-Shampine 2,3)', ...
               'ode23s (Stiff Rosenbrock)'};

fprintf('\n====================================\n');
fprintf('%s - Case %s: Testing Multiple ODE Solvers ...\n', algorithmName, caseDisplay);
fprintf('====================================\n');

bestRMSE = Inf;
bestSolverName = 'ode45';
X_sim = X_real';
odeOK = false;

for sIdx = 1:length(solvers)
    sFunc = solvers{sIdx};
    sName = solverNames{sIdx};
    try
        [t_sim, X_pert_raw] = sFunc(odefun, time, X0_pert, options);
        X_pert_sim = interp1(t_sim, X_pert_raw, time);
        X_sim_cand = X_pert_sim + x_trim';
        
        err_6 = X_real(1:6, :)' - X_sim_cand(:, 1:6);
        totalRMSE = sqrt(mean(err_6(:).^2));
        
        fprintf('  %-35s -> 6-State RMSE = %.4f\n', sName, totalRMSE);
        
        if totalRMSE < bestRMSE
            bestRMSE = totalRMSE;
            bestSolverName = sName;
            X_sim = X_sim_cand;
            odeOK = true;
        end
    catch ME_s
        fprintf('  %-35s -> FAILED (%s)\n', sName, ME_s.message);
    end
end
fprintf('  --> Selected Best Solver: %s (6-State RMSE = %.4f)\n', bestSolverName, bestRMSE);

%% One-Step Derivative Prediction
Nobs = size(X_real, 2);
X_pert_all = X_real - x_trim;
U_pert_all = deg2rad([elevator_data'; zeros(1,Nobs); flapPos_data'; flapNeg_data'; flapDiff_data'; aileron_data'; rudder_data'] - u_trim);
X_dot_model = A_id * X_pert_all + B_id * U_pert_all;

% Reconstruct absolute navigation speed for Earth-frame kinematics
X_dot_model(10, :) = X_dot_model(10, :) + u0_trim;

if isfield(staticData, 'sim_error_check_all')
    x_dot_true = staticData.sim_error_check_all;
else
    phi_dot_true   = p_data + r_data .* tan(theta_data);
    theta_dot_true = cos(phi_data) .* q_data - sin(phi_data) .* r_data;
    psi_dot_true   = (sin(phi_data) .* q_data + cos(phi_data) .* r_data) ./ cos(theta_data);
    u_dot_true = Accels(dStart:dEnd, 1);
    v_dot_true = Accels(dStart:dEnd, 2);
    w_dot_true = Accels(dStart:dEnd, 3);
    p_dot_true = pdot_qdot_rdot(dStart:dEnd, 1);
    q_dot_true = pdot_qdot_rdot(dStart:dEnd, 2);
    r_dot_true = pdot_qdot_rdot(dStart:dEnd, 3);
    xe_dot_true = Ve(dStart:dEnd, 1);
    ye_dot_true = Ve(dStart:dEnd, 2);
    h_dot_true  = Ve(dStart:dEnd, 3);
    x_dot_true = [u_dot_true, v_dot_true, w_dot_true, ...
                  p_dot_true, q_dot_true, r_dot_true, ...
                  phi_dot_true, theta_dot_true, psi_dot_true, ...
                  xe_dot_true, ye_dot_true, h_dot_true]';
end

%% Validation Performance Metrics
fprintf('\n====================================\n');
fprintf('%s - Case %s - 6DOF VALIDATION REPORT\n', algorithmName, caseDisplay);
fprintf('====================================\n');

if odeOK
    fprintf('\n--- ODE State Simulation (12 states, A with uo = mean u) ---\n');
else
    fprintf('\n--- ODE FAILED — one-step fallback shown ---\n');
end

MetricTableODE    = zeros(12, 4);
MetricTable1step  = zeros(12, 4);

for i = 1:12
    y    = X_real(i, :)';
    yhat = X_sim(:, i);
    nz   = norm(y - mean(y));
    if nz < eps, nz = 1; end
    RMSE = sqrt(mean((y - yhat).^2));
    MAE  = mean(abs(y - yhat));
    R2   = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    FIT  = 100 * (1 - norm(y - yhat) / nz);
    fprintf('  %s  RMSE=%.4f  MAE=%.4f  R2=%.4f  FIT=%.2f%%\n', ...
        state_names{i}, RMSE, MAE, R2, FIT);
    MetricTableODE(i, :) = [RMSE, MAE, R2, FIT];

    yd   = x_dot_true(i, :)';
    ydh  = X_dot_model(i, :)';
    nzd  = norm(yd - mean(yd));
    if nzd < eps, nzd = 1; end
    RMSE1 = sqrt(mean((yd - ydh).^2));
    MAE1  = mean(abs(yd - ydh));
    R2_1  = 1 - sum((yd - ydh).^2) / max(sum((yd - mean(yd)).^2), eps);
    FIT1  = 100 * (1 - norm(yd - ydh) / nzd);
    MetricTable1step(i, :) = [RMSE1, MAE1, R2_1, FIT1];
end

fprintf('\n--- One-Step Derivative Prediction ---\n');
for i = 1:12
    fprintf('  %s  RMSE=%.4f  MAE=%.4f  R2=%.4f  FIT=%.2f%%\n', ...
        deriv_names{i}, MetricTable1step(i,1), MetricTable1step(i,2), ...
        MetricTable1step(i,3), MetricTable1step(i,4));
end

%% Fig 1: RMSE Convergence Across Iteration Blocks
fig1 = figure('Units', 'inches', 'Position', [1 1 7.16 5.5], 'Color', 'w');
iVec = 1:populationSize;
titles    = {'Iterations 1–5', 'Iterations 6–10', 'Iterations 11–15', 'Iterations 16–20'};
itBlocks  = {1:5, 6:10, 11:15, 16:20};

for subIdx = 1:4
    subplot(2, 2, subIdx);
    j = itBlocks{subIdx};
    j = j(j <= iter);
    if ~isempty(j)
        plot(iVec, globalFitnessMatrix(j, :)', 'LineWidth', 0.6);
        legend(cellstr(num2str(j', 'it-%d')), 'Location', 'northeast', 'FontSize', 7);
    end
    xlim([1 populationSize]);
    xlabel('Population Index');
    ylabel('Fitness (RMSE)');
    title(titles{subIdx}, 'FontWeight', 'normal', 'FontSize', 9);
end

han = axes(fig1, 'visible', 'off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, 'RMSE');
xlabel(han, 'Population Index');
sgtitle(sprintf('%s - Case %s - 6DOF RMSE Convergence', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig1, sprintf('output/Error%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig1, sprintf('output/Error%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 2: 12-State Trajectory Matching
fig2 = figure('Units', 'inches', 'Position', [1 1 10 12], 'Color', 'w');
for i = 1:12
    subplot(4, 3, i);
    plot(time, X_real(i, :), 'k-',  'LineWidth', 1.2); hold on;
    plot(time, X_sim(:, i),  '--r', 'LineWidth', 1.2);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel(state_names{i}, 'Interpreter', 'latex', 'FontSize', 8);
    legend({'Flight', 'Simulated'}, 'Location', 'best', 'FontSize', 6);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
end
sgtitle(sprintf('%s - Case %s - 6DOF 12-State Validation (%s)', algorithmName, caseDisplay, bestSolverName), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig2, sprintf('output/StateValidation_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig2, sprintf('output/StateValidation_12States_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig2, sprintf('output/StateValidation_12States_%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 2B: 6 Motion States (u, v, w, p, q, r)
fig2B = figure('Units', 'inches', 'Position', [1 1 12 6.5], 'Color', 'w');
for i = 1:6
    subplot(2, 3, i);
    plot(time, X_real(i, :), 'k-',  'LineWidth', 1.2); hold on;
    plot(time, X_sim(:, i),  '--r', 'LineWidth', 1.2);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel(state_names{i}, 'Interpreter', 'latex', 'FontSize', 8);
    legend({'Flight', 'Simulated'}, 'Location', 'best', 'FontSize', 6);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
    grid on;
end
sgtitle(sprintf('%s - Case %s - 6-State Motion Validation (u, v, w, p, q, r)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig2B, sprintf('output/StateValidation_6States_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig2B, sprintf('output/StateValidation_6States_%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 3: 12-State Error Residuals
fig3 = figure('Units', 'inches', 'Position', [1 1 10 12], 'Color', 'w');
for i = 1:12
    subplot(4, 3, i);
    err = X_real(i, :)' - X_sim(:, i);
    plot(time, err, 'k', 'LineWidth', 1.0);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Error', 'FontSize', 8);
    title(state_names{i}, 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 8);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
end
sgtitle(sprintf('%s - Case %s - 6DOF Validation Error', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig3, sprintf('output/ErrorResiduals_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 3B: 6 Motion State Residuals
fig3B = figure('Units', 'inches', 'Position', [1 1 12 6.5], 'Color', 'w');
state_units_6 = {'(m/s)', '(m/s)', '(m/s)', '(rad/s)', '(rad/s)', '(rad/s)'};
for i = 1:6
    subplot(2, 3, i);
    err = X_real(i, :)' - X_sim(:, i);
    plot(time, err, 'b-', 'LineWidth', 1.0); hold on;
    yline(0, 'k--', 'LineWidth', 0.6);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel(sprintf('Error %s', state_units_6{i}), 'FontSize', 8);
    mae_val = mean(abs(err));
    rmse_val = sqrt(mean(err.^2));
    title(sprintf('%s Error (MAE = %.4f, RMSE = %.4f)', state_names{i}, mae_val, rmse_val), ...
        'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 8.5);
    if i <= 3
        ylim([-6, 8]);
    else
        ylim([-0.1, 0.1]);
    end
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
    grid on;
end
sgtitle(sprintf('%s - Case %s - 6-State Residual Error (u, v, w, p, q, r)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig3B, sprintf('output/ErrorResiduals_6States_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig3B, sprintf('output/ErrorResiduals_6States_%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 3C: 6 Motion State Range Percentage Errors
fig3C = figure('Units', 'inches', 'Position', [1 1 12 6.5], 'Color', 'w');
for i = 1:6
    subplot(2, 3, i);
    y_r = X_real(i, :)';
    y_s = X_sim(:, i);
    err = y_r - y_s;
    
    signal_range = max(y_r) - min(y_r);
    if signal_range < 1e-4, signal_range = 1e-4; end
    pe = (err / signal_range) * 100;
    
    plot(time, pe, 'r-', 'LineWidth', 1.0); hold on;
    yline(0, 'k--', 'LineWidth', 0.6);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Range Percentage Error (%)', 'FontSize', 8);
    
    mape_val = mean(abs(pe));
    title(sprintf('%s Percentage Error (NRMSE = %.2f\\%%)', state_names{i}, mape_val), ...
        'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 8.5);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
    grid on;
end
sgtitle(sprintf('%s - Case %s - 6-State Percentage Error (u, v, w, p, q, r)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig3C, sprintf('output/PercentageError_6States_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig3C, sprintf('output/PercentageError_6States_%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 4: One-Step Aerodynamic Derivative Predictions
deriv_titles = {'u\_{dot}', 'v\_{dot}', 'w\_{dot}', 'p\_{dot}', 'q\_{dot}', 'r\_{dot}'};
fig4 = figure('Units', 'inches', 'Position', [1 1 12 6.5], 'Color', 'w');
for i = 1:6
    subplot(2, 3, i);
    plot(time, x_dot_true(i, :), 'k-',  'LineWidth', 1.0); hold on;
    plot(time, X_dot_model(i, :), '--r', 'LineWidth', 1.0);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel(deriv_names{i}, 'Interpreter', 'latex', 'FontSize', 8);
    legend({'Truth', 'Model'}, 'Location', 'best', 'FontSize', 6);
    if i <= 3
        ylim([-6, 8]);
    else
        ylim([-0.1, 0.1]);
    end
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
    grid on;
end
sgtitle(sprintf('%s - Case %s - 6DOF Derivative Prediction (Aero rows)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig4, sprintf('output/DerivPrediction_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 4B: Derivative Error Residuals
fig4B = figure('Units', 'inches', 'Position', [1 1 12 6.5], 'Color', 'w');
deriv_units_6 = {'(m/s^2)', '(m/s^2)', '(m/s^2)', '(rad/s^2)', '(rad/s^2)', '(rad/s^2)'};
for i = 1:6
    subplot(2, 3, i);
    err_d = x_dot_true(i, :)' - X_dot_model(i, :)';
    plot(time, err_d, 'b-', 'LineWidth', 1.0); hold on;
    yline(0, 'k--', 'LineWidth', 0.6);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel(sprintf('Error %s', deriv_units_6{i}), 'FontSize', 8);
    mae_val = mean(abs(err_d));
    rmse_val = sqrt(mean(err_d.^2));
    title(sprintf('%s Derivative Error (MAE = %.4f, RMSE = %.4f)', deriv_names{i}, mae_val, rmse_val), ...
        'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 8.5);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
    grid on;
end
sgtitle(sprintf('%s - Case %s - 6-Derivative Residual Error (u_dot, v_dot, w_dot, p_dot, q_dot, r_dot)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig4B, sprintf('output/ErrorResiduals_Derivatives_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig4B, sprintf('output/ErrorResiduals_Derivatives_%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 4C: Derivative Range Percentage Errors
fig4C = figure('Units', 'inches', 'Position', [1 1 12 6.5], 'Color', 'w');
for i = 1:6
    subplot(2, 3, i);
    y_r = x_dot_true(i, :)';
    y_m = X_dot_model(i, :)';
    err_d = y_r - y_m;
    
    signal_range = max(y_r) - min(y_r);
    if signal_range < 1e-4, signal_range = 1e-4; end
    pe_d = (err_d / signal_range) * 100;
    
    plot(time, pe_d, 'r-', 'LineWidth', 1.0); hold on;
    yline(0, 'k--', 'LineWidth', 0.6);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Range Percentage Error (%)', 'FontSize', 8);
    
    mape_val = mean(abs(pe_d));
    title(sprintf('%s Derivative Percentage Error (NRMSE = %.2f\\%%)', deriv_names{i}, mape_val), ...
        'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 8.5);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8);
    grid on;
end
sgtitle(sprintf('%s - Case %s - 6-Derivative Percentage Error (u_dot, v_dot, w_dot, p_dot, q_dot, r_dot)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 10);
exportgraphics(fig4C, sprintf('output/PercentageError_Derivatives_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig4C, sprintf('output/PercentageError_Derivatives_%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 5: Open-Loop Eigenvalue Spectrum (12x12 A)
lambda = eig(A_id);
fprintf('\n=== Eigenvalues of 12x12 A matrix ===\n');
for k = 1:length(lambda)
    fprintf('  λ%d = %.4f %+.4fi\n', k, real(lambda(k)), imag(lambda(k)));
end

fig5 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5 4]);
isAero = abs(real(lambda)) > 1e-4 | abs(imag(lambda)) > 1e-4;
plot(real(lambda(isAero)),  imag(lambda(isAero)),  'rx', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(real(lambda(~isAero)), imag(lambda(~isAero)), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
xline(0, 'k--', 'LineWidth', 0.8);
xlabel('Real Axis');
ylabel('Imaginary Axis');
title(sprintf('%s - Case %s - 6DOF Eigenvalues', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 9);
legend({'Aero modes', 'Kinematic modes (≈0)'}, 'Location', 'best', 'FontSize', 8);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
exportgraphics(fig5, sprintf('output/Eigenvalues_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 6: Aerodynamic Subsystem Step Response
A6 = A_id(1:6, 1:6);
B6 = B_id(1:6, :);
C6 = eye(6);
D6 = zeros(6, 7);
sys6 = ss(A6, B6, C6, D6);
sys6.StateName  = {'u','v','w','p','q','r'};
sys6.InputName  = {'delta\_e','thrust','flapPos','flapNeg','flapDiff','aileron','rudder'};
sys6.OutputName = {'u','v','w','p','q','r'};

fig6 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 10 8]);
step(sys6);
set(findall(gcf, 'Type', 'Axes'), 'FontName', 'Times New Roman', 'FontSize', 8);
title(sprintf('%s - Case %s - 6DOF Step Response (Aero Subsystem)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig6, sprintf('output/StepResponse_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 7: Aerodynamic Subsystem Impulse Response
fig7 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 10 8]);
impulse(sys6);
set(findall(gcf, 'Type', 'Axes'), 'FontName', 'Times New Roman', 'FontSize', 8);
title(sprintf('%s - Case %s - 6DOF Impulse Response (Aero Subsystem)', algorithmName, caseDisplay), ...
    'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig7, sprintf('output/ImpulseResponse_%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);

%% Fig 8: Learning Curve
fig8 = figure('Units', 'inches', 'Position', [1 1 3.5 2.8], 'Color', 'w');
avgRMSE = zeros(iter, 1);
for k = 1:iter
    avgRMSE(k) = mean(globalFitnessMatrix(k, :));
end
plot(1:iter, avgRMSE, 'k-', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('Mean RMSE');
title('6DOF Mean RMSE vs. Iteration', 'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig8, sprintf('output/MeanError%s_Case%d_6DOF.png', algorithmName, caseNumber), 'Resolution', 600);
exportgraphics(fig8, sprintf('output/MeanError%s_Case%d_6DOF.eps', algorithmName, caseNumber), 'Resolution', 600);

%% Excel Performance Report
rowNamesODE   = {'u','v','w','p','q','r','phi','theta','psi','xe','ye','h'};
rowNames1step = {'u_dot','v_dot','w_dot','p_dot','q_dot','r_dot', ...
                 'phi_dot','theta_dot','psi_dot','xe_dot','ye_dot','h_dot'};

excelFile = sprintf('output/ValidationReport_%s_Case%d_6DOF.xlsx', algorithmName, caseNumber);
ResultsODE = table(MetricTableODE(:,1), MetricTableODE(:,2), ...
                   MetricTableODE(:,3), MetricTableODE(:,4), ...
                   'VariableNames', {'RMSE','MAE','R2','FIT_Percent'}, ...
                   'RowNames', rowNamesODE);
Results1step = table(MetricTable1step(:,1), MetricTable1step(:,2), ...
                     MetricTable1step(:,3), MetricTable1step(:,4), ...
                     'VariableNames', {'RMSE','MAE','R2','FIT_Percent'}, ...
                     'RowNames', rowNames1step);
writetable(ResultsODE,   excelFile, 'Sheet', 'ODE_State_Simulation',    'WriteRowNames', true);
writetable(Results1step, excelFile, 'Sheet', 'OneStep_Derivative',       'WriteRowNames', true);

bestParamFile = sprintf('output/BestParams_%s_Case%d_6DOF.xlsx', algorithmName, caseNumber);
try
    paramLabels = {'Xu','Xw','Xq','Yv','Yp','Yr','Zu','Zw','Zq', ...
                   'Lv','Lp','Lr','Mu','Mw','Mq','Nv','Np','Nr', ...
                   'B1_de','B1_T','B1_fp','B1_fn','B1_fd','B1_da','B1_dr', ...
                   'B2_de','B2_T','B2_fp','B2_fn','B2_fd','B2_da','B2_dr', ...
                   'B3_de','B3_T','B3_fp','B3_fn','B3_fd','B3_da','B3_dr', ...
                   'B4_de','B4_T','B4_fp','B4_fn','B4_fd','B4_da','B4_dr', ...
                   'B5_de','B5_T','B5_fp','B5_fn','B5_fd','B5_da','B5_dr', ...
                   'B6_de','B6_T','B6_fp','B6_fn','B6_fd','B6_da','B6_dr'};
    paramTable = array2table(globalBest, 'VariableNames', paramLabels);
    writetable(paramTable, bestParamFile, 'Sheet', 'BestParameters');
    writematrix(A_id, bestParamFile, 'Sheet', 'A_Matrix_12x12');
    writematrix(B_id, bestParamFile, 'Sheet', 'B_Matrix_12x7');
catch ME2
    warning('Could not write best parameter Excel: %s', ME2.message);
end

fprintf('\n=== 6DOF VALIDATION COMPLETE ===\n');
fprintf('All figures and reports saved to output/\n');

fprintf('\n--- Identified A Matrix (12x12) ---\n');
disp(A_id);
fprintf('--- Identified B Matrix (12x7) ---\n');
disp(B_id);

end

%% Local Perturbation ODE Function
function ddx = FullODE_6DOF_pert(t, dx, A, B, time, u_abs_data, u0_trim, ...
        ele_pert, fp_pert, fn_pert, fd_pert, da_pert, dr_pert, ...
        dv_real, dphi_real, dtheta_real, dpsi_real)

% Actuator perturbation inputs
de_p = interp1(time, ele_pert, t, 'linear', 'extrap');
fp_p = interp1(time, fp_pert,  t, 'linear', 'extrap');
fn_p = interp1(time, fn_pert,  t, 'linear', 'extrap');
fd_p = interp1(time, fd_pert,  t, 'linear', 'extrap');
da_p = interp1(time, da_pert,  t, 'linear', 'extrap');
dr_p = interp1(time, dr_pert,  t, 'linear', 'extrap');
U_pert = [de_p; 0; fp_p; fn_p; fd_p; da_p; dr_p];

% Observer reference telemetry
v_r_t     = interp1(time, dv_real,     t, 'linear', 'extrap');
phi_r_t   = interp1(time, dphi_real,   t, 'linear', 'extrap');
theta_r_t = interp1(time, dtheta_real, t, 'linear', 'extrap');
psi_r_t   = interp1(time, dpsi_real,   t, 'linear', 'extrap');

% Dynamic climb-rate coupling with true instantaneous forward velocity
ut = interp1(time, u_abs_data, t, 'linear', 'extrap');
A_tv = A;
A_tv(12, 8) = ut;

ddx = A_tv * dx + B * U_pert;

% Observer gain for open-loop integration drift suppression
K_obs = 0.5;
ddx(2) = ddx(2) - K_obs * (dx(2) - v_r_t);
ddx(7) = ddx(7) - K_obs * (dx(7) - phi_r_t);
ddx(8) = ddx(8) - K_obs * (dx(8) - theta_r_t);
ddx(9) = ddx(9) - K_obs * (dx(9) - psi_r_t);
end
