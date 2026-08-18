function ALO24(inputData)
% ALO24  Ant Lion Optimizer for 24-parameter longitudinal aircraft parameter estimation.
% Identifies stability and multi-surface control derivatives (elevators, multi-segment flaps).

arguments
    inputData InputData
end

%% Flight data and problem configuration
staticData = PrepareFlightData24(inputData);

% Unpack kinematics, accelerations, and surface deflections for post-processing
Vb = staticData.Vb;
pqr = staticData.pqr;
phi_theta_psi = staticData.phi_theta_psi;
Xe = staticData.Xe;
Accels = staticData.Accels;
elevator = staticData.elevator;
flapPos = staticData.flapPos;
flapNeg = staticData.flapNeg;
flapDiff = staticData.flapDiff;
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
globalrmse = [];

%% Initialize population and best solution tracking
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
globalBest = population(1,:);
globalBestFitness = inf;

%% Main optimization loop
for iter = 1:maxIterations

    % Cost function evaluation across the swarm
    fitnessValues = EvaluateFitness24(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx,:);

    % Update global best solution
    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, :) = fitnessValues;
    
    % Early convergence check on mean cost plateau
    globalrmse(iter) = mean(fitnessValues);
    if iter > 1
        diff_val = abs(globalrmse(iter) - globalrmse(iter-1));
        if diff_val < 0.0001
            disp("algorithm end condition reached");
            break;
        end
    end
    
    % Position update bounded by parameter limits
    population = UpdatePositions(population, globalBest, lowerBounds, upperBounds);

end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('Execution time: %.5f seconds\n', elapsedTime);

%% Extract identified state-space model and export data
u = staticData.u_all(end);
[A, B] = formatParameters24(bestSolution, u);
disp('Best solution found. Matrices');
disp("A Matrix")
disp(A);
disp("B Matrix")
disp(B);

save("ALO-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
