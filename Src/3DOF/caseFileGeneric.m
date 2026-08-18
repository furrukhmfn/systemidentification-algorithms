function caseFileGeneric(caseNumber, populationSize, offSetIncreament, iteration, algorithmFunc, algorithmName)
% Generic runner for 24-parameter longitudinal model identification using any algorithm handle.

if nargin < 6 || isempty(algorithmName)
    s = functions(algorithmFunc);
    algorithmName = s.function;
end

% Add algorithm library path if needed
algoDir = fullfile(fileparts(mfilename('fullpath')), 'more_algorithms');
if exist(algoDir, 'dir')
    addpath(algoDir);
end

if ~exist('output', 'dir')
    mkdir('output');
end

% Telemetry files and output paths
accelerationFileName = "cruise_acceleration.mat";
outputStatesFileName = "cruise_outputStates.mat";
accutatorsFilesName = "cruise_acctuators.mat";
workSpaceFileName = fullfile('output', algorithmName + "_Case" + num2str(caseNumber) + "_24.mat");
loggingFileName = fullfile('output', algorithmName + "_Logging_Case" + num2str(caseNumber) + "_24.txt");

parameterCount = 24;

% Stability and control parameter bounds
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

% Execute optimization routine
feval(algorithmFunc, inputData);

diary off;

end
