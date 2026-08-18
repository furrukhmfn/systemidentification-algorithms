function ALO6DOF(inputData)
% ALO6DOF  Ant Lion Optimizer for 6-DOF aircraft parameter identification.
% Estimates 60 stability and control derivatives for a 12-state linear model:
%   States (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
%   Inputs  (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]

arguments
    inputData InputData
end

% Load and pre-process flight telemetry
staticData = PrepareFlightData6DOF(inputData);

Vb             = staticData.Vb;
pqr            = staticData.pqr;
phi_theta_psi  = staticData.phi_theta_psi;
Xe             = staticData.Xe;
Accels         = staticData.Accels;
aileron        = staticData.aileron;
rudder         = staticData.rudder;
elevator       = staticData.elevator;
flapPos        = staticData.flapPos;
flapNeg        = staticData.flapNeg;
flapDiff       = staticData.flapDiff;
Ve             = staticData.Ve;
pdot_qdot_rdot = staticData.pdot_qdot_rdot;

lowerBounds   = inputData.LowerBounds;
upperBounds   = inputData.UpperBounds;
maxIterations = inputData.MaxIteration;
populationSize = inputData.PopulationSize;
parameterCount = inputData.ParameterCount;
valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);

% Initialize population within parameter bounds
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);

globalBest = population(1, :);
globalBestFitness = inf;

%% Optimization Loop
for iter = 1:maxIterations

    fitnessValues = EvaluateFitness6DOF(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx, :);

    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, :) = fitnessValues;

    % Update positions toward elite solution
    population = UpdatePositions(population, globalBest, lowerBounds, upperBounds);
end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('ALO6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== ALO6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('ALO-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
