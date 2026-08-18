function staticData = PrepareFlightData(inputData)
% Loads flight telemetry logs and extracts longitudinal state/derivative 
% arrays for model identification and fitness evaluation.

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
elevator = accutatorData.accutators.Data(offset:limit,2);
Ve = outputData.output_states.Ve.Data(offset:limit,:);

% Vehicle mass and nominal thrust
Kg2lb = 2.20462;
massBody = 22932 / Kg2lb;
massFuel = 2948 / Kg2lb;
mass = massBody + massFuel;
thrust = mass * 4;

% Longitudinal states: [u, w, q, theta, h]
u_all = Vb(:, 1);
w_all = Vb(:, 3);
q_all = pqr(:, 2);
theta_all = phi_theta_psi(2, :)';
h_all = -Xe(:, 3);
elevator_all = elevator(:, 1);

% State derivative references from accelerometers and rate gyros
u_dot_all = Accels(:, 1);
w_dot_all = Accels(:, 3);
q_dot_all = pdot_qdot_rdot(:, 2);
h_dot_all = Ve(:, 3);
phi_all = phi_theta_psi(1, :)';

% Kinematic pitch rate transformation: theta_dot = q*cos(phi) - r*sin(phi)
theta_dot_all = cos(phi_all) .* q_all - sin(phi_all) .* pqr(:, 3);
sim_error_check_all = [u_dot_all, w_dot_all, q_dot_all, theta_dot_all, h_dot_all]';

% Package precomputed arrays and objective weights
staticData.u_all = u_all;
staticData.w_all = w_all;
staticData.q_all = q_all;
staticData.theta_all = theta_all;
staticData.h_all = h_all;
staticData.elevator_all = elevator_all;
staticData.thrust = 0;
staticData.sim_error_check_all = sim_error_check_all;
staticData.Weights = [0.5; 0.2; 0.1; 0.1; 0.1];
staticData.valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

% Telemetry signals retained for validation and plotting
staticData.Vb = Vb;
staticData.pqr = pqr;
staticData.phi_theta_psi = phi_theta_psi;
staticData.Xe = Xe;
staticData.Accels = Accels;
staticData.elevator = elevator;
staticData.Ve = Ve;
staticData.pdot_qdot_rdot = pdot_qdot_rdot;

end
