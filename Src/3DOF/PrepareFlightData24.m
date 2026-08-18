function staticData = PrepareFlightData24(inputData)
% Loads flight telemetry logs and extracts longitudinal state, control, and derivative 
% arrays (including flap channels) for 24-parameter system identification.

% Load telemetry logs
accelData = load(inputData.AccelerationFileName);
outputData = load(inputData.OutputStatesFilesName);
accutatorData = load(inputData.AccutatorsFileName);

populationNumber = inputData.PopulationNumber;
offset = inputData.OffSet;
populationSize = inputData.PopulationSize;
limit = populationNumber + offset;

% Slice flight telemetry for the designated sample window
pdot_qdot_rdot = outputData.output_states.pdot_qdot_rdot.Data(offset:limit,:);
pqr = outputData.output_states.pqr.Data(offset:limit,:);

phi_theta_psi = squeeze(outputData.output_states.phi_theta_psi.Data);
phi_theta_psi = phi_theta_psi(:, offset:limit);

Vb = outputData.output_states.Vb.Data(offset:limit,:);
Xe = outputData.output_states.Xe.Data(offset:limit,:);
Accels = accelData.acceleration_bb.Data(offset:limit,:);
Ve = outputData.output_states.Ve.Data(offset:limit,:);

% Actuator channels (elevator and asymmetric/symmetric flaps)
elevator = accutatorData.accutators.Data(offset:limit,2);
flapPos = accutatorData.accutators.Data(offset:limit,4);
flapNeg = accutatorData.accutators.Data(offset:limit,5);
flapDiff = accutatorData.accutators.Data(offset:limit,6);

% Vehicle mass properties
Kg2lb = 2.20462;
massBody = 22932 / Kg2lb;
massFuel = 2948 / Kg2lb;
mass = massBody + massFuel;
thrust = mass * 4;

% Trim transient boundary margins
margin = inputData.ValidationMargin;
idxS = 1 + margin;
idxE = size(Vb, 1) - margin;

% Longitudinal states: [u, w, q, theta, h]
u_all = Vb(idxS:idxE, 1);
w_all = Vb(idxS:idxE, 3);
q_all = pqr(idxS:idxE, 2);
theta_all = phi_theta_psi(2, idxS:idxE)';
h_all = -Xe(idxS:idxE, 3);

% Primary and secondary control surface deflections
elevator_all = elevator(idxS:idxE, 1);
flapPos_all = flapPos(idxS:idxE, 1);
flapNeg_all = flapNeg(idxS:idxE, 1);
flapDiff_all = flapDiff(idxS:idxE, 1);

% Reference acceleration and derivative signals
u_dot_all = Accels(idxS:idxE, 1);
w_dot_all = Accels(idxS:idxE, 3);
q_dot_all = pdot_qdot_rdot(idxS:idxE, 2);
h_dot_all = Ve(idxS:idxE, 3);
phi_all = phi_theta_psi(1, idxS:idxE)';

% Kinematic pitch rate transformation: theta_dot = q*cos(phi) - r*sin(phi)
theta_dot_all = cos(phi_all) .* q_all - sin(phi_all) .* pqr(idxS:idxE, 3);
sim_error_check_all = [u_dot_all, w_dot_all, q_dot_all, theta_dot_all, h_dot_all]';
time_all = outputData.output_states.pdot_qdot_rdot.Time(offset:limit);
time_all = time_all(idxS:idxE) - time_all(idxS);

% Package precomputed arrays and objective weights
staticData.u_all = u_all;
staticData.w_all = w_all;
staticData.q_all = q_all;
staticData.theta_all = theta_all;
staticData.h_all = h_all;
staticData.elevator_all = elevator_all;
staticData.flapPos_all = flapPos_all;
staticData.flapNeg_all = flapNeg_all;
staticData.flapDiff_all = flapDiff_all;
staticData.thrust = 0;
staticData.time_all = time_all;
staticData.sim_error_check_all = sim_error_check_all;
staticData.Weights = [0.5, 0.2, 0.2 0.05 0.05]';
fprintf('  Weights: [%.4f, %.4f, %.4f, %.4f, %.4f]\n', staticData.Weights);
staticData.valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

% Telemetry signals retained for validation and plotting
staticData.Vb = Vb;
staticData.pqr = pqr;
staticData.phi_theta_psi = phi_theta_psi;
staticData.Xe = Xe;
staticData.Accels = Accels;
staticData.elevator = elevator;
staticData.flapPos = flapPos;
staticData.flapNeg = flapNeg;
staticData.flapDiff = flapDiff;
staticData.Ve = Ve;
staticData.pdot_qdot_rdot = pdot_qdot_rdot;

end
