function JADE(inputData)
% JADE  Adaptive metaheuristic optimization for 24-parameter longitudinal model identification.
% Estimates aircraft stability and multi-control derivatives from flight test data.

arguments
    inputData InputData
end

%% Flight data and problem setup
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

b = 1;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

%% Initialize search population
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
bestFlames = population;

%% Main optimization loop
for iter = 1:maxIterations
    % Cost function evaluation
    fitnessValues = EvaluateFitness24(population, staticData);

    [sortedFitnessValues, sortedIndices] = sort(fitnessValues);
    bestFlames = population(sortedIndices, :);

    globalFitnessMatrix(iter, :) = fitnessValues;

    % Early convergence check
    globalrmse(iter) = mean(fitnessValues);
    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("algorithm end condition reached");
        break;
    end

    % Spiral trajectory position update toward elite candidates
    t = (iter / maxIterations) * 2 * pi;
    newPopulation = zeros(populationSize, parameterCount);
    for i = 1:populationSize
        for j = 1:parameterCount
            D = abs(bestFlames(i, j) - population(i, j));
            newPopulation(i, j) = D * exp(b * t) * cos(t) + bestFlames(i, j);
        end
        newPopulation(i, :) = max(lowerBounds, min(upperBounds, newPopulation(i, :)));
    end
    population = newPopulation;
end

globalBest = bestFlames(1, :);
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

save("JADE-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
