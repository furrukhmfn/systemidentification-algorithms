function caseFile(caseNumber, populationSize, offSetIncreament, iteration)
% Runs Ant Lion Optimizer (ALO) for the 15-parameter longitudinal flight model.

if ~exist('output', 'dir')
    mkdir('output');
end

% Telemetry files and output paths
accelerationFileName = "cruise_acceleration.mat";
outputStatesFileName = "cruise_outputStates.mat";
accutatorsFilesName = "cruise_acctuators.mat";
workSpaceFileName = fullfile('output', "ALO_Case"+num2str(caseNumber)+".mat");
loggingFileName = fullfile('output', "ALO_Logging_Case"+num2str(caseNumber)+".txt");

% 15-parameter longitudinal stability and control bounds
parameterCount = 15;
lowerBounds = [-10,0,0,-10,0,0,0,-20,-20,0,0,0,0,-20,-20];
upperBounds = [0,100,5,0,10,10,0,20,20,10,10,10,10,10,10];

maxIterations = iteration;

diary(loggingFileName);

inputData = InputData(accelerationFileName, outputStatesFileName, accutatorsFilesName, workSpaceFileName,...
    lowerBounds, upperBounds, maxIterations, populationSize, parameterCount);
if(offSetIncreament ~= 0)
    inputData.IncreamentOffSet(offSetIncreament);
end

% Execute optimization
ALO(inputData);

% Post-processing and response graphics
inputGraphicsForLongDynmaics = InputGraphicsData(workSpaceFileName, "ALO", caseNumber);
CreateGraphicsLongitudinalDynamics(inputGraphicsForLongDynmaics);

diary off;

end