function [fitnessValues, x_dot] = EvaluateFitness(population, staticData)
% Evaluates candidate parameter vectors (15-parameter longitudinal model) 
% using weighted RMSE against measured state derivatives.

% Unpack dimensional stability and control derivatives
p1 = population(:, 1);
p2 = population(:, 2);
p3 = population(:, 3);
p4 = population(:, 4);
p5 = population(:, 5);
p6 = population(:, 6);
p7 = population(:, 7);
p8 = population(:, 8);
p9 = population(:, 9);
p10 = population(:, 10);
p11 = population(:, 11);
p12 = population(:, 12);
p13 = population(:, 13);
p14 = population(:, 14);
p15 = population(:, 15);

% Flight states and trim constants
u_all = staticData.u_all;
w_all = staticData.w_all;
q_all = staticData.q_all;
theta_all = staticData.theta_all;
elevator_all = staticData.elevator_all;
thrust = staticData.thrust;
g = staticData.valueOfGravitationConstant;

M = size(population, 1); % population size
N = size(u_all, 1);       % flight data sample count

% Regressor matrix [u, w, q, delta_e, thrust]
R = [u_all, w_all, q_all, elevator_all, thrust * ones(N, 1)];

% Longitudinal force and moment derivative blocks
Q1 = [p1, p2, p3, p10, p11]; % u-dot equation parameters
Q2 = [p4, p5, p6, p12, p13]; % w-dot equation parameters
Q3 = [p7, p8, p9, p14, p15]; % q-dot equation parameters

% Predicted state derivatives across trajectory [M x N]
x_dot_1 = Q1 * R' - g * theta_all';                     % u_dot
x_dot_2 = Q2 * R';                                      % w_dot
x_dot_3 = Q3 * R';                                      % q_dot
x_dot_4 = ones(M, 1) * q_all';                          % theta_dot (kinematic)
x_dot_5 = ones(M, 1) * (-w_all + u_all .* theta_all)';  % h_dot (kinematic)

% Flight test reference derivatives
u_dot_true = staticData.sim_error_check_all(1, :);
w_dot_true = staticData.sim_error_check_all(2, :);
q_dot_true = staticData.sim_error_check_all(3, :);
theta_dot_true = staticData.sim_error_check_all(4, :);
h_dot_true = staticData.sim_error_check_all(5, :);

% State MSE per candidate across trajectory
mse_1 = mean((x_dot_1 - u_dot_true).^2, 2);
mse_2 = mean((x_dot_2 - w_dot_true).^2, 2);
mse_3 = mean((x_dot_3 - q_dot_true).^2, 2);
mse_4 = mean((x_dot_4 - theta_dot_true).^2, 2);
mse_5 = mean((x_dot_5 - h_dot_true).^2, 2);

mse_matrix = [mse_1, mse_2, mse_3, mse_4, mse_5];

% Weighted trajectory RMSE
fitnessValues = sqrt(mse_matrix * staticData.Weights);

% Diagnostic derivative snapshot [5 x M] for backwards compatibility
idx = mod((0:M-1)', N) + 1;
u_all_trunc = u_all(idx);
w_all_trunc = w_all(idx);
q_all_trunc = q_all(idx);
theta_all_trunc = theta_all(idx);
elevator_all_trunc = elevator_all(idx);

x_dot_1_diag = p1 .* u_all_trunc + p2 .* w_all_trunc + p3 .* q_all_trunc - g .* theta_all_trunc + p10 .* elevator_all_trunc + p11 .* thrust;
x_dot_2_diag = p4 .* u_all_trunc + p5 .* w_all_trunc + p6 .* q_all_trunc + p12 .* elevator_all_trunc + p13 .* thrust;
x_dot_3_diag = p7 .* u_all_trunc + p8 .* w_all_trunc + p9 .* q_all_trunc + p14 .* elevator_all_trunc + p15 .* thrust;
x_dot_4_diag = q_all_trunc;
x_dot_5_diag = -w_all_trunc + u_all_trunc .* theta_all_trunc;

x_dot = [x_dot_1_diag, x_dot_2_diag, x_dot_3_diag, x_dot_4_diag, x_dot_5_diag]';

end
