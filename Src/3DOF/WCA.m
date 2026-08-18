function WCA(inputData)
% WCA  Water Cycle Algorithm for 24-parameter longitudinal aircraft parameter estimation.
% Models river and stream flows toward the sea (global best) with evaporation and rainfall dispersion.

arguments
    inputData InputData
end

%% Flight data and optimization setup
staticData = PrepareFlightData24(inputData);

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

nRivers = 4;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

%% Initialize search population
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
globalBest = population(1, :);
globalBestFitness = inf;

%% Main optimization loop
for iter = 1:maxIterations
    % Evaluate candidate model fitness
    fitnessValues = EvaluateFitness24(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx, :);
    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, :) = fitnessValues;

    % Early convergence check
    globalrmse(iter) = mean(fitnessValues);
    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("algorithm end condition reached");
        break;
    end

    % Rank population to designate sea and rivers
    [~, sortIndex] = sort(fitnessValues);
    population = population(sortIndex, :);

    sea = population(1, :);
    rivers = population(2:nRivers + 1, :);

    % Stream flow towards rivers and sea with evaporation/rainfall effect
    for i = nRivers + 2:populationSize
        riverIdx = randi([1 nRivers]);
        flowDistance = 0.05 * rand() .* (rivers(riverIdx, :) - population(i, :));
        population(i, :) = population(i, :) + flowDistance;
        if rand() < 0.1
            population(riverIdx, :) = population(riverIdx, :) + 0.01 * rand() .* (sea - rivers(riverIdx, :));
        end
        population(i, :) = max(lowerBounds, min(upperBounds, population(i, :)));
    end
end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('Execution time: %.5f seconds\n', elapsedTime);

%% Extract state-space matrices and export
u = staticData.u_all(end);
[A, B] = formatParameters24(bestSolution, u);
disp('Best solution found. Matrices');
disp("A Matrix")
disp(A);
disp("B Matrix")
disp(B);

save("WCA-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
