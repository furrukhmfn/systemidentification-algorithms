function SSA(inputData)
% SSA  Salp Swarm Algorithm for 24-parameter longitudinal aircraft parameter estimation.
% Identifies stability and control derivatives via leader-follower chain propagation.

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

    % Exploration-exploitation balance parameter
    c1 = 2 * exp(-((4 * iter / maxIterations)^2));
    newPopulation = zeros(populationSize, parameterCount);
    for i = 1:populationSize
        if i == 1
            % Leader salp position update towards global best
            for j = 1:parameterCount
                c2 = rand(); c3 = rand();
                if c3 < 0.5
                    newPopulation(i, j) = globalBest(j) + c1 * ((upperBounds(j) - lowerBounds(j)) * c2 + lowerBounds(j));
                else
                    newPopulation(i, j) = globalBest(j) - c1 * ((upperBounds(j) - lowerBounds(j)) * c2 + lowerBounds(j));
                end
            end
        else
            % Follower salp chain update
            newPopulation(i, :) = (population(i, :) + population(i - 1, :)) / 2;
        end
        newPopulation(i, :) = max(lowerBounds, min(upperBounds, newPopulation(i, :)));
    end
    population = newPopulation;
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

save("SSA-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
