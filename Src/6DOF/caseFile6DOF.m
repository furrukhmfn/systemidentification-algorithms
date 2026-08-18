function caseFile6DOF(caseNumber, populationSize, offSetIncreament, iteration, varargin)
% caseFile6DOF  Sets up and runs 6-DOF aircraft parameter identification via GWO.
%
% Identifies a 12-state, 7-input state-space model:
%   States (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
%   Inputs  (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]
%   Params (60): 18 stability derivatives + 42 control derivatives

% Add algorithm library to search path
algoDir = fullfile(fileparts(mfilename('fullpath')), 'more_algorithms');
if exist(algoDir, 'dir')
    addpath(algoDir);
end

if ~exist('output', 'dir')
    mkdir('output');
end

%% Telemetry Configuration
% Default identification and out-of-sample validation partitions
train_accelFile    = 'acceleration.mat';
train_statesFile   = 'outputStates.mat';
train_actuatorFile = 'acctuators.mat';
train_offset       = 1;
train_length       = 4000;

test_accelFile     = 'cruise_acceleration.mat';
test_statesFile    = 'cruise_outputStates.mat';
test_actuatorFile  = 'cruise_acctuators.mat';
test_offset        = 4001;
test_length        = 4000;

% Optional override via custom telemetry structs
if nargin >= 5 && isstruct(varargin{1})
    train_accelFile    = varargin{1}.accel;
    train_statesFile   = varargin{1}.states;
    train_actuatorFile = varargin{1}.actuators;
    if isfield(varargin{1}, 'offset'), train_offset = varargin{1}.offset; end
    if isfield(varargin{1}, 'length'), train_length = varargin{1}.length; end
end
if nargin >= 6 && isstruct(varargin{2})
    test_accelFile     = varargin{2}.accel;
    test_statesFile    = varargin{2}.states;
    test_actuatorFile  = varargin{2}.actuators;
    if isfield(varargin{2}, 'offset'), test_offset = varargin{2}.offset; end
    if isfield(varargin{2}, 'length'), test_length = varargin{2}.length; end
end

workSpaceFileName = fullfile('output', sprintf('GWO_Case%d_6DOF.mat', caseNumber));
loggingFileName   = fullfile('output', sprintf('GWO_Logging_Case%d_6DOF.txt', caseNumber));

parameterCount = 60;

%% Parameter Search Bounds
% 18 stability derivatives + 42 control derivatives
% Physical constraints:
%   Xu <= 0  (speed damping / drag)
%   Yv <= 0  (sideforce opposing sideslip)
%   Zu <= 0  (lift increase with axial speed)
%   Zw < 0   (heave damping)
%   Lp <= 0  (roll damping)
%   Mw <= 0  (static longitudinal stability)
%   Mq <= 0  (pitch damping)
%   Nr <= 0  (yaw damping)
stabLower = [-10, -100, -10,  -10, -10, -10,  -10, -10, -10,  -20, -20, -10,  -1, -20, -20,    0, -10, -20];
stabUpper = [  0,  100,  10, -0.15,  10,  10,    0, -eps, 10,    0, -eps, 10,   1,   0, -eps,  20,  10, -eps];

% Control derivatives (6 aerodynamic rows x 7 inputs)
ctrlLower = -20 * ones(1, 42);
ctrlUpper =  20 * ones(1, 42);

lowerBounds = [stabLower, ctrlLower];
upperBounds = [stabUpper, ctrlUpper];

%% Run Identification
diary(loggingFileName);

inputData = InputData(train_accelFile, train_statesFile, train_actuatorFile, ...
    workSpaceFileName, lowerBounds, upperBounds, iteration, populationSize, parameterCount);

inputData.OffSet = train_offset;
inputData.PopulationNumber = train_length;

% Attach validation telemetry
inputData = inputData.SetTestData(test_accelFile, test_statesFile, test_actuatorFile, test_offset, test_length);

if offSetIncreament ~= 0
    inputData = inputData.IncreamentOffSet(offSetIncreament);
end

fprintf('caseFile6DOF: Running GWO6DOF — Case %d, popSize=%d, iter=%d\n', ...
    caseNumber, populationSize, iteration);

GWO6DOF(inputData);

% Generate validation figures and metrics
CreateGraphics6DOF(workSpaceFileName, caseNumber, 'GWO');

diary off;

end
