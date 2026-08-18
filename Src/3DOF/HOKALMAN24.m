function HOKALMAN24(inputData)
% HOKALMAN24  Subspace state-space identification (N4SID / Ho-Kalman) for longitudinal aircraft dynamics.
% Identifies an unconstrained state-space model from I/O flight data and projects the simulated
% trajectories onto the structured 24-parameter longitudinal equations of motion via least squares.

arguments
    inputData InputData
end

%% Flight data and parameters
staticData = PrepareFlightData24(inputData);

Vb = staticData.Vb;
pqr = staticData.pqr;
phi_theta_psi = staticData.phi_theta_psi;
Xe = staticData.Xe;
elevator = staticData.elevator;
flapPos = staticData.flapPos;
flapNeg = staticData.flapNeg;
flapDiff = staticData.flapDiff;

lowerBounds = inputData.LowerBounds;
upperBounds = inputData.UpperBounds;
valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

tic;
fprintf("HOKALMAN24: Subspace identification (N4SID)\n");

u_all = staticData.u_all;
w_all = staticData.w_all;
q_all = staticData.q_all;
theta_all = staticData.theta_all;
elevator_all = staticData.elevator_all;
flapPos_all = staticData.flapPos_all;
flapNeg_all = staticData.flapNeg_all;
flapDiff_all = staticData.flapDiff_all;
N = size(u_all, 1);
thrust = staticData.thrust;
g = valueOfGravitationConstant;
dt = mean(diff(staticData.time_all));

% Assemble measured output and input matrices
Y = [u_all, w_all, q_all, theta_all, -Xe(:, 3)];
U_mat = [elevator_all, thrust * ones(N, 1), flapPos_all, flapNeg_all, flapDiff_all];

%% Step 1: Subspace state-space identification (N4SID)
dataObj = iddata(Y, U_mat, dt);
opt = n4sidOptions('Focus', 'simulation', 'N4Weight', 'CVA', 'Display', 'off');
sys = n4sid(dataObj, 5, opt);

%% Step 2: Linear simulation of identified subspace model
Y_sim = lsim(sys, U_mat, staticData.time_all, Y(1, :));

%% Step 3: Compute state time derivatives from simulated states
u_sim = Y_sim(:, 1);
w_sim = Y_sim(:, 2);
q_sim = Y_sim(:, 3);
theta_sim = Y_sim(:, 4);
h_sim = Y_sim(:, 5);

du_sim = gradient(u_sim, dt);
dw_sim = gradient(w_sim, dt);
dq_sim = gradient(q_sim, dt);
dtheta_sim = q_sim;                      % Pitch rate kinematic relation
dh_sim = -w_sim + u_sim .* theta_sim;    % Climb rate kinematic relation

%% Step 4: Project onto structured 24-parameter aircraft model via least squares
% Regressor matrix containing longitudinal states and control inputs
R = [u_sim, w_sim, q_sim, elevator_all, thrust * ones(N, 1), ...
     flapPos_all, flapNeg_all, flapDiff_all];

theta1 = R \ (du_sim + g * theta_sim);
theta2 = R \ dw_sim;
theta3 = R \ dq_sim;

p = zeros(1, 24);
p(1:9) = [theta1(1:3)', theta2(1:3)', theta3(1:3)'];
p(10:15) = [theta1(4:5)', theta2(4:5)', theta3(4:5)'];
p(16:24) = [theta1(6:8)', theta2(6:8)', theta3(6:8)'];
p = max(lowerBounds, min(upperBounds, p));

bestSolution = p;
globalBest = p;
uo = u_all(end);
[A, B] = formatParameters24(p, uo);

elapsedTime = toc;
fprintf('Execution time: %.5f seconds\n', elapsedTime);

disp('Identified 24 parameters (projected from N4SID):');
fprintf('  Xu=%.4f  Xw=%.4f  Xq=%.4f\n', p(1), p(2), p(3));
fprintf('  Zu=%.4f  Zw=%.4f  Zq=%.4f\n', p(4), p(5), p(6));
fprintf('  Mu=%.4f  Mw=%.4f  Mq=%.4f\n', p(7), p(8), p(9));
disp("A Matrix")
disp(A);
disp("B Matrix")
disp(B);

% Evaluate cost of the structured projected model
fit_n4 = EvaluateFitness24(p, staticData);
fprintf('Structured model RMSE from N4SID projection: %.4f\n', fit_n4);

% Variables for post-processing and plotting compatibility
fitnessValues = fit_n4;
globalrmse = fitnessValues;
globalBestFitness = fitnessValues;
iter = 1;
populationSize = 1;
globalFitnessMatrix = fitnessValues;

save("HOKALMAN-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
