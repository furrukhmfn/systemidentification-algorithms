function staticData = PrepareFlightData6DOF(inputData)
% PrepareFlightData6DOF  Extracts and pre-processes 6-DOF flight telemetry for identification.
% Formulates state and control perturbation variables around cruise trim:
%   States (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
%   Inputs  (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]

% Load raw telemetry MAT files
accelData     = load(inputData.AccelerationFileName);
outputData    = load(inputData.OutputStatesFilesName);
accutatorData = load(inputData.AccutatorsFileName);

populationNumber = inputData.PopulationNumber;
offset = inputData.OffSet;
limit  = populationNumber + offset;

% Verify telemetry length bounds
dataSize = size(outputData.output_states.pdot_qdot_rdot.Data, 1);
if limit > dataSize
    warning('PrepareFlightData6DOF: limit=%d exceeds data length=%d — clamping to %d.', ...
        limit, dataSize, dataSize);
    limit = dataSize;
end

% Extract slice intervals
pdot_qdot_rdot = outputData.output_states.pdot_qdot_rdot.Data(offset:limit, :);
pqr            = outputData.output_states.pqr.Data(offset:limit, :);

phi_theta_psi  = squeeze(outputData.output_states.phi_theta_psi.Data);
phi_theta_psi  = phi_theta_psi(:, offset:limit);

Vb     = outputData.output_states.Vb.Data(offset:limit, :);
Xe     = outputData.output_states.Xe.Data(offset:limit, :);
Accels = accelData.acceleration_bb.Data(offset:limit, :);
Ve     = outputData.output_states.Ve.Data(offset:limit, :);

% Actuator channels: 1=aileron, 2=rudder, 3=elevator, 4=flapPos, 5=flapNeg, 6=flapDiff
aileron  = accutatorData.accutators.Data(offset:limit, 1);
rudder   = accutatorData.accutators.Data(offset:limit, 2);
elevator = accutatorData.accutators.Data(offset:limit, 3);
flapPos  = accutatorData.accutators.Data(offset:limit, 4);
flapNeg  = accutatorData.accutators.Data(offset:limit, 5);
flapDiff = accutatorData.accutators.Data(offset:limit, 6);

% Windowing with validation margin
margin = inputData.ValidationMargin;
idxS = 1 + margin;
idxE = size(Vb, 1) - margin;

%% Flight State Variables
u_all     = Vb(idxS:idxE, 1);               % body-axis axial velocity (m/s)
v_all     = Vb(idxS:idxE, 2);               % body-axis lateral velocity (m/s)
w_all     = Vb(idxS:idxE, 3);               % body-axis normal velocity (m/s)
p_all     = pqr(idxS:idxE, 1);              % roll rate (rad/s)
q_all     = pqr(idxS:idxE, 2);              % pitch rate (rad/s)
r_all     = pqr(idxS:idxE, 3);              % yaw rate (rad/s)
phi_all   = phi_theta_psi(1, idxS:idxE)';   % bank angle (rad)
theta_all = phi_theta_psi(2, idxS:idxE)';   % pitch angle (rad)
psi_all   = phi_theta_psi(3, idxS:idxE)';   % heading angle (rad)
xe_all    = Xe(idxS:idxE, 1);               % North position (m)
ye_all    = Xe(idxS:idxE, 2);               % East position (m)
h_all     = -Xe(idxS:idxE, 3);              % Altitude above MSL (m)

%% Trim Equilibrium States
u0_t   = mean(u_all);    v0_t = mean(v_all);    w0_t = mean(w_all);
p0_t   = mean(p_all);    q0_t = mean(q_all);    r0_t = mean(r_all);
phi0_t = mean(phi_all);  theta0_t = mean(theta_all);  psi0_t = mean(psi_all);
xe0_t  = mean(xe_all);   ye0_t = mean(ye_all);  h0_t = mean(h_all);

fprintf('  Trim: u0=%.2f v0=%.4f w0=%.4f p0=%.6f q0=%.6f r0=%.6f theta0=%.4f\n', ...
    u0_t, v0_t, w0_t, p0_t, q0_t, r0_t, theta0_t);

% State perturbations (Delta_x = x - x_trim)
du_all     = u_all - u0_t;
dv_all     = v_all - v0_t;
dw_all     = w_all - w0_t;
dp_all     = p_all - p0_t;
dq_all     = q_all - q0_t;
dr_all     = r_all - r0_t;
dphi_all   = phi_all   - phi0_t;
dtheta_all = theta_all - theta0_t;

%% Reference Flight State Derivatives
% Aerodynamic acceleration derivatives from body accelerometers and gyro rates
u_dot_all = Accels(idxS:idxE, 1);           % du/dt (m/s^2)
v_dot_all = Accels(idxS:idxE, 2);           % dv/dt (m/s^2)
w_dot_all = Accels(idxS:idxE, 3);           % dw/dt (m/s^2)
p_dot_all = pdot_qdot_rdot(idxS:idxE, 1);   % dp/dt (rad/s^2)
q_dot_all = pdot_qdot_rdot(idxS:idxE, 2);   % dq/dt (rad/s^2)
r_dot_all = pdot_qdot_rdot(idxS:idxE, 3);   % dr/dt (rad/s^2)

% Euler angle kinematic rates
phi_dot_all   = p_all + r_all .* tan(theta_all);
theta_dot_all = cos(phi_all) .* q_all - sin(phi_all) .* r_all;
psi_dot_all   = (sin(phi_all) .* q_all + cos(phi_all) .* r_all) ./ cos(theta_all);

% Earth-frame navigation velocity rates
xe_dot_all = Ve(idxS:idxE, 1);              % dxe/dt (m/s)
ye_dot_all = Ve(idxS:idxE, 2);              % dye/dt (m/s)
h_dot_all  = Ve(idxS:idxE, 3);              % dh/dt  (m/s)

%% Actuator Telemetry and Perturbations
aileron_all  = aileron(idxS:idxE, 1);
rudder_all   = rudder(idxS:idxE, 1);
elevator_all = elevator(idxS:idxE, 1);
flapPos_all  = flapPos(idxS:idxE, 1);
flapNeg_all  = flapNeg(idxS:idxE, 1);
flapDiff_all = flapDiff(idxS:idxE, 1);

da0_t  = mean(aileron_all);
dr0_t  = mean(rudder_all);
de0_t  = mean(elevator_all);
fp0_t  = mean(flapPos_all);
fn0_t  = mean(flapNeg_all);
fd0_t  = mean(flapDiff_all);

% Actuator deflections converted to radians for standard SI derivatives
dda_all  = deg2rad(aileron_all  - da0_t);
ddr_all  = deg2rad(rudder_all   - dr0_t);
dde_all  = deg2rad(elevator_all - de0_t);
dfp_all  = deg2rad(flapPos_all  - fp0_t);
dfn_all  = deg2rad(flapNeg_all  - fn0_t);
dfd_all  = deg2rad(flapDiff_all - fd0_t);

x_trim = [u0_t; v0_t; w0_t; p0_t; q0_t; r0_t; ...
          phi0_t; theta0_t; psi0_t; xe0_t; ye0_t; h0_t];
u_trim = [de0_t; 0; fp0_t; fn0_t; fd0_t; da0_t; dr0_t];

fprintf('  Actuator trim: da0=%.4f dr0=%.4f de0=%.4f fp0=%.4f fn0=%.4f fd0=%.4f\n', ...
    da0_t, dr0_t, de0_t, fp0_t, fn0_t, fd0_t);

%% Normalization Factors
% Standard deviation floor to avoid division by zero on low-energy channels
AeroDerivStd = max([std(u_dot_all), std(v_dot_all), std(w_dot_all), ...
                    std(p_dot_all), std(q_dot_all), std(r_dot_all)]', 1e-4);
fprintf('  6DOF aero deriv std: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', AeroDerivStd);

sim_error_check_all = [u_dot_all, v_dot_all, w_dot_all, ...
                        p_dot_all, q_dot_all, r_dot_all, ...
                        phi_dot_all, theta_dot_all, psi_dot_all, ...
                        xe_dot_all,  ye_dot_all,  h_dot_all]';

time_all = outputData.output_states.pdot_qdot_rdot.Time(offset:limit);
time_all = time_all(idxS:idxE) - time_all(idxS);

% Relative weighting for 6 aerodynamic states (kinematic rows weighted 0)
Weights = [0.20; 0.20; 0.20; 0.15; 0.15; 0.10; 0; 0; 0; 0; 0; 0];
fprintf('  6DOF NMSE weights: [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f] | kinematic: 0\n', ...
    Weights(1), Weights(2), Weights(3), Weights(4), Weights(5), Weights(6));

%% Output Data Struct
staticData.u_all     = u_all;
staticData.v_all     = v_all;
staticData.w_all     = w_all;
staticData.p_all     = p_all;
staticData.q_all     = q_all;
staticData.r_all     = r_all;
staticData.phi_all   = phi_all;
staticData.theta_all = theta_all;
staticData.psi_all   = psi_all;
staticData.xe_all    = xe_all;
staticData.ye_all    = ye_all;
staticData.h_all     = h_all;

staticData.aileron_all  = aileron_all;
staticData.rudder_all   = rudder_all;
staticData.elevator_all = elevator_all;
staticData.flapPos_all  = flapPos_all;
staticData.flapNeg_all  = flapNeg_all;
staticData.flapDiff_all = flapDiff_all;
staticData.thrust = 0;

staticData.sim_error_check_all = sim_error_check_all;
staticData.time_all = time_all;
staticData.Weights  = Weights;
staticData.AeroDerivStd = AeroDerivStd;
staticData.valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

staticData.du_all     = du_all;
staticData.dv_all     = dv_all;
staticData.dw_all     = dw_all;
staticData.dp_all     = dp_all;
staticData.dq_all     = dq_all;
staticData.dr_all     = dr_all;
staticData.dphi_all   = dphi_all;
staticData.dtheta_all = dtheta_all;
staticData.dda_all    = dda_all;
staticData.ddr_all    = ddr_all;
staticData.dde_all    = dde_all;
staticData.dfp_all    = dfp_all;
staticData.dfn_all    = dfn_all;
staticData.dfd_all    = dfd_all;

staticData.x_trim  = x_trim;
staticData.u_trim  = u_trim;
staticData.u0_trim = u0_t;

% Telemetry matrices retained for plotting and validation
staticData.Vb             = Vb;
staticData.pqr            = pqr;
staticData.phi_theta_psi  = phi_theta_psi;
staticData.Xe             = Xe;
staticData.Accels         = Accels;
staticData.aileron        = aileron;
staticData.rudder         = rudder;
staticData.elevator       = elevator;
staticData.flapPos        = flapPos;
staticData.flapNeg        = flapNeg;
staticData.flapDiff       = flapDiff;
staticData.Ve             = Ve;
staticData.pdot_qdot_rdot = pdot_qdot_rdot;

end
