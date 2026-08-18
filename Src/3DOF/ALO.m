function ALO(inputData)
% ALO  Ant Lion Optimizer for 3-DOF longitudinal aircraft system identification.
% Estimates stability and control derivatives from flight test time-history data.

arguments
    inputData InputData
end

%% Flight data and problem setup
staticData = PrepareFlightData(inputData);

% Unpack state and control time-histories for workspace export
Vb = staticData.Vb;
pqr = staticData.pqr;
phi_theta_psi = staticData.phi_theta_psi;
Xe = staticData.Xe;
Accels = staticData.Accels;
elevator = staticData.elevator;
Ve = staticData.Ve;
pdot_qdot_rdot = staticData.pdot_qdot_rdot;

lowerBounds = inputData.LowerBounds;
upperBounds = inputData.UpperBounds;
maxIterations = inputData.MaxIteration;
populationSize = inputData.PopulationSize;
parameterCount = inputData.ParameterCount;
valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);

%% Initialize population and best solution tracking
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
globalBest = population(1,:);
globalBestFitness = inf;

%% Main optimization loop
for iter = 1:maxIterations

    % Cost function evaluation across the swarm
    fitnessValues = EvaluateFitness(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx,:);

    % Update global best candidate
    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, :) = fitnessValues;

    % Position update bounded by search space limits
    population = UpdatePositions(population, globalBest, lowerBounds, upperBounds);

end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('Execution time: %.5f seconds\n', elapsedTime);

%% Format identified state-space matrices and save results
u = staticData.u_all(end);
[A,B] = formatParameters(bestSolution, u, valueOfGravitationConstant);
disp('Best solution found. Matrices');
disp("A Matrix")
disp(A);
disp("B Matrix")
disp(B);
[dirPath, fileName, fileExt] = fileparts(inputData.WorkspaceFileName);
save(fullfile(dirPath, "BasicVariables_" + fileName + fileExt), "globalFitnessMatrix", "globalBest", "valueOfGravitationConstant", "fitnessValues", "elapsedTime");
save(inputData.WorkspaceFileName);

end