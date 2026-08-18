function caseFile24(caseNumber, populationSize, offSetIncreament, iteration)
% Runs Ant Lion Optimizer (ALO24) for the 24-parameter longitudinal flight model.

if ~exist('output', 'dir')
    mkdir('output');
end

% Telemetry files and output paths
accelerationFileName = "cruise_acceleration.mat";
outputStatesFileName = "cruise_outputStates.mat";
accutatorsFilesName = "cruise_acctuators.mat";
workSpaceFileName = fullfile('output', "ALO_Case"+num2str(caseNumber)+"_24.mat");
loggingFileName = fullfile('output', "ALO_Logging_Case"+num2str(caseNumber)+"_24.txt");

parameterCount = 24;

% Search space bounds for 24 parameters:
% Stability derivatives: [Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq]
% Primary controls:      [X_de, X_dT, Z_de, Z_dT, M_de, M_dT]
% Flap controls:         [X_fp, X_fn, X_fd, Z_fp, Z_fn, Z_fd, M_fp, M_fn, M_fd]
% Physical constraints: Zw < 0 (heave damping), Mw <= 0 (static stability), Mq <= 0 (pitch damping)
lowerBounds = [-10, -100, -10, -10, -10, -10, -1, -20, -20, ...
    -20, -20, -20, -20, -20, -20, ...
    -20, -20, -20, -20, -20, -20, -20, -20, -20];
upperBounds = [ 0,  100,  10,   0,   -eps, 10,  1,   0,   0, ...
     20,  20,  20,  20,  20,  20, ...
     20,  20,  20,  20,  20,  20,  20,  20,  20];

maxIterations = iteration;

diary(loggingFileName);

inputData = InputData(accelerationFileName, outputStatesFileName, accutatorsFilesName, workSpaceFileName,...
    lowerBounds, upperBounds, maxIterations, populationSize, parameterCount);

if offSetIncreament ~= 0
    inputData.IncreamentOffSet(offSetIncreament);
end

% Execute optimization
ALO24(inputData);

% Post-processing and response validation graphics
CreateGraphicsLongitudinalDynamics24(workSpaceFileName, caseNumber, "ALO");

diary off;

end
