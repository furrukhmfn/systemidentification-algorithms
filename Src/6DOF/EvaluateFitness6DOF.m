function [fitnessValues, x_dot, rmse_per_state] = EvaluateFitness6DOF(population, staticData)
% EvaluateFitness6DOF  Calculates variance-normalized NMSE cost for 6-DOF identification.
% Evaluates candidate 60-parameter vectors using linearized perturbation dynamics:
%   Delta_x_dot = A * Delta_x + B * Delta_u
%
% States (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
% Inputs  (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]
% Params (60): 18 stability derivatives + 42 control derivatives

% State and actuator perturbations around trim
du_all     = staticData.du_all;
dv_all     = staticData.dv_all;
dw_all     = staticData.dw_all;
dp_all     = staticData.dp_all;
dq_all     = staticData.dq_all;
dr_all     = staticData.dr_all;
dphi_all   = staticData.dphi_all;
dtheta_all = staticData.dtheta_all;
dda_all    = staticData.dda_all;
ddr_all    = staticData.ddr_all;
dde_all    = staticData.dde_all;
dfp_all    = staticData.dfp_all;
dfn_all    = staticData.dfn_all;
dfd_all    = staticData.dfd_all;

% Absolute states for kinematic expressions
p_all     = staticData.p_all;
q_all     = staticData.q_all;
r_all     = staticData.r_all;
phi_all   = staticData.phi_all;
theta_all = staticData.theta_all;
u_all     = staticData.u_all;
w_all     = staticData.w_all;

g = staticData.valueOfGravitationConstant;
W = staticData.Weights;
uo = staticData.u0_trim;

M = size(population, 1);
N = size(du_all, 1);

T_col = zeros(N, 1);

%% Regressors in Perturbation Coordinates (N x 10)
% Inputs: [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]
R_long = [du_all, dw_all, dq_all, dde_all, T_col, dfp_all, dfn_all, dfd_all, dda_all, ddr_all];
R_lat  = [dv_all, dp_all, dr_all, dde_all, T_col, dfp_all, dfn_all, dfd_all, dda_all, ddr_all];

%% Parameter Subsets (M x 10)
Q1 = population(:, [1,  2,  3,  19, 20, 21, 22, 23, 24, 25]);  % u_dot
Q2 = population(:, [4,  5,  6,  26, 27, 28, 29, 30, 31, 32]);  % v_dot
Q3 = population(:, [7,  8,  9,  33, 34, 35, 36, 37, 38, 39]);  % w_dot
Q4 = population(:, [10, 11, 12, 40, 41, 42, 43, 44, 45, 46]);  % p_dot
Q5 = population(:, [13, 14, 15, 47, 48, 49, 50, 51, 52, 53]);  % q_dot
Q6 = population(:, [16, 17, 18, 54, 55, 56, 57, 58, 59, 60]);  % r_dot

%% Aerodynamic State Derivatives (M x N)
xd_u = Q1 * R_long';
xd_v = Q2 * R_lat';
xd_w = Q3 * R_long';
xd_p = Q4 * R_lat';
xd_q = Q5 * R_long';
xd_r = Q6 * R_lat';

%% Kinematic Rates (analytical formulas)
phi_dot_kin   = ones(M,1) * (p_all + r_all .* tan(theta_all))';
theta_dot_kin = ones(M,1) * (cos(phi_all) .* q_all - sin(phi_all) .* r_all)';
psi_dot_kin   = ones(M,1) * ((sin(phi_all) .* q_all + cos(phi_all) .* r_all) ./ cos(theta_all))';
xed_kin       = ones(M,1) * staticData.sim_error_check_all(10, :);
yed_kin       = ones(M,1) * staticData.sim_error_check_all(11, :);
hd_kin        = ones(M,1) * (-w_all + u_all .* theta_all)';

%% Reference Flight Derivatives (12 x N)
u_dot_true     = staticData.sim_error_check_all(1,  :);
v_dot_true     = staticData.sim_error_check_all(2,  :);
w_dot_true     = staticData.sim_error_check_all(3,  :);
p_dot_true     = staticData.sim_error_check_all(4,  :);
q_dot_true     = staticData.sim_error_check_all(5,  :);
r_dot_true     = staticData.sim_error_check_all(6,  :);
phi_dot_true   = staticData.sim_error_check_all(7,  :);
theta_dot_true = staticData.sim_error_check_all(8,  :);
psi_dot_true   = staticData.sim_error_check_all(9,  :);
xed_true       = staticData.sim_error_check_all(10, :);
yed_true       = staticData.sim_error_check_all(11, :);
hd_true        = staticData.sim_error_check_all(12, :);

%% Mean Squared Errors per State
mse1  = mean((xd_u         - u_dot_true    ).^2, 2);
mse2  = mean((xd_v         - v_dot_true    ).^2, 2);
mse3  = mean((xd_w         - w_dot_true    ).^2, 2);
mse4  = mean((xd_p         - p_dot_true    ).^2, 2);
mse5  = mean((xd_q         - q_dot_true    ).^2, 2);
mse6  = mean((xd_r         - r_dot_true    ).^2, 2);
mse7  = mean((phi_dot_kin  - phi_dot_true  ).^2, 2);
mse8  = mean((theta_dot_kin- theta_dot_true).^2, 2);
mse9  = mean((psi_dot_kin  - psi_dot_true  ).^2, 2);
mse10 = mean((xed_kin      - xed_true      ).^2, 2);
mse11 = mean((yed_kin      - yed_true      ).^2, 2);
mse12 = mean((hd_kin       - hd_true       ).^2, 2);

mse_matrix = [mse1, mse2, mse3, mse4, mse5, mse6, ...
              mse7, mse8, mse9, mse10, mse11, mse12];

%% Normalized Cost Function (NMSE)
% Balances translational (m/s^2) and angular (rad/s^2) acceleration residual scales
ADS2 = (staticData.AeroDerivStd).^2;

nmse_aero = [mse1./ADS2(1), mse2./ADS2(2), mse3./ADS2(3), ...
             mse4./ADS2(4), mse5./ADS2(5), mse6./ADS2(6)];

W6 = staticData.Weights(1:6);
fitnessValues = sqrt(nmse_aero * W6);

rmse_per_state = sqrt(mse_matrix);

%% Open-Loop Stability Regularization
% Penalize unstable real eigenvalue components
penalty = zeros(M, 1);
for i = 1:M
    pr = population(i, :);
    A_i = [pr(1),  0,     pr(2),  0,     pr(3),      0,     0,   -g,  0, 0, 0, 0;
            0,    pr(4),   0,    pr(5),   0,   pr(6)-uo,    g,    0,  0, 0, 0, 0;
           pr(7),  0,     pr(8),  0,   pr(9)+uo,     0,     0,    0,  0, 0, 0, 0;
            0,   pr(10),  0,   pr(11),   0,    pr(12),   0,    0,  0, 0, 0, 0;
          pr(13),  0,   pr(14),  0,    pr(15),       0,     0,    0,  0, 0, 0, 0;
            0,   pr(16),  0,   pr(17),   0,    pr(18),   0,    0,  0, 0, 0, 0;
            0,    0,      0,    1,       0,          0,     0,    0,  0, 0, 0, 0;
            0,    0,      0,    0,       1,          0,     0,    0,  0, 0, 0, 0;
            0,    0,      0,    0,       0,          1,     0,    0,  0, 0, 0, 0;
            1,    0,      0,    0,       0,          0,     0,    0,  0, 0, 0, 0;
            0,    1,      0,    0,       0,          0,     0,    0,  0, 0, 0, 0;
            0,    0,     -1,    0,       0,          0,     0,   uo,  0, 0, 0, 0];

    lam = eig(A_i);
    maxReal = max(real(lam));
    if maxReal > 0.001
        penalty(i) = 10 + 100 * maxReal;
    end
end
fitnessValues = fitnessValues + penalty;

%% One-Step Diagnostic Predictions
idx = mod((0:M-1)', N) + 1;

du_t   = du_all(idx);     dv_t   = dv_all(idx);     dw_t   = dw_all(idx);
dp_t   = dp_all(idx);     dq_t   = dq_all(idx);     dr_t   = dr_all(idx);
dtheta_t = dtheta_all(idx);  dphi_t = dphi_all(idx);
dda_t  = dda_all(idx);
ddr_t  = ddr_all(idx);
dde_t  = dde_all(idx);
dfp_t  = dfp_all(idx);    dfn_t  = dfn_all(idx);    dfd_t  = dfd_all(idx);

p_t    = p_all(idx);   q_t    = q_all(idx);   r_t    = r_all(idx);
phi_t  = phi_all(idx); theta_t = theta_all(idx);
u_t    = u_all(idx);   w_t    = w_all(idx);
T_t    = zeros(M, 1);

p1  = population(:,1);  p2  = population(:,2);  p3  = population(:,3);
p4  = population(:,4);  p5  = population(:,5);  p6  = population(:,6);
p7  = population(:,7);  p8  = population(:,8);  p9  = population(:,9);
p10 = population(:,10); p11 = population(:,11); p12 = population(:,12);
p13 = population(:,13); p14 = population(:,14); p15 = population(:,15);
p16 = population(:,16); p17 = population(:,17); p18 = population(:,18);
p19 = population(:,19); p20 = population(:,20); p21 = population(:,21);
p22 = population(:,22); p23 = population(:,23); p24 = population(:,24);
p25 = population(:,25);
p26 = population(:,26); p27 = population(:,27); p28 = population(:,28);
p29 = population(:,29); p30 = population(:,30); p31 = population(:,31);
p32 = population(:,32);
p33 = population(:,33); p34 = population(:,34); p35 = population(:,35);
p36 = population(:,36); p37 = population(:,37); p38 = population(:,38);
p39 = population(:,39);
p40 = population(:,40); p41 = population(:,41); p42 = population(:,42);
p43 = population(:,43); p44 = population(:,44); p45 = population(:,45);
p46 = population(:,46);
p47 = population(:,47); p48 = population(:,48); p49 = population(:,49);
p50 = population(:,50); p51 = population(:,51); p52 = population(:,52);
p53 = population(:,53);
p54 = population(:,54); p55 = population(:,55); p56 = population(:,56);
p57 = population(:,57); p58 = population(:,58); p59 = population(:,59);
p60 = population(:,60);

xd1  = p1.*du_t  + p2.*dw_t  + p3.*dq_t  - g.*dtheta_t + p19.*dde_t + p20.*T_t + p21.*dfp_t + p22.*dfn_t + p23.*dfd_t + p24.*dda_t + p25.*ddr_t;
xd2  = p4.*dv_t  + p5.*dp_t  + p6.*dr_t  + g.*dphi_t   + p26.*dde_t + p27.*T_t + p28.*dfp_t + p29.*dfn_t + p30.*dfd_t + p31.*dda_t + p32.*ddr_t;
xd3  = p7.*du_t  + p8.*dw_t  + p9.*dq_t               + p33.*dde_t + p34.*T_t + p35.*dfp_t + p36.*dfn_t + p37.*dfd_t + p38.*dda_t + p39.*ddr_t;
xd4  = p10.*dv_t + p11.*dp_t + p12.*dr_t               + p40.*dde_t + p41.*T_t + p42.*dfp_t + p43.*dfn_t + p44.*dfd_t + p45.*dda_t + p46.*ddr_t;
xd5  = p13.*du_t + p14.*dw_t + p15.*dq_t               + p47.*dde_t + p48.*T_t + p49.*dfp_t + p50.*dfn_t + p51.*dfd_t + p52.*dda_t + p53.*ddr_t;
xd6  = p16.*dv_t + p17.*dp_t + p18.*dr_t               + p54.*dde_t + p55.*T_t + p56.*dfp_t + p57.*dfn_t + p58.*dfd_t + p59.*dda_t + p60.*ddr_t;

xd7  = p_t;
xd8  = cos(phi_t) .* q_t - sin(phi_t) .* r_t;
xd9  = (sin(phi_t) .* q_t + cos(phi_t) .* r_t) ./ cos(theta_t);
xd10 = u_t;
xd11 = q_t;
xd12 = -w_t + u_t .* theta_t;

x_dot = [xd1, xd2, xd3, xd4, xd5, xd6, xd7, xd8, xd9, xd10, xd11, xd12]';

end
