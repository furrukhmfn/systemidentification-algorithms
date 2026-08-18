function WOA(inputData)
% WOA  Whale Optimization Algorithm for 24-parameter longitudinal aircraft parameter estimation.
% Estimates stability and control derivatives via encircling, bubble-net spiral, and search mechanisms.

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

    % Linearly decrease exploration coefficient
    a = 2 - iter * (2 / maxIterations);
    newPopulation = zeros(populationSize, parameterCount);
    for i = 1:populationSize
        r1 = rand(); r2 = rand();
        A = 2 * a * r1 - a;
        C = 2 * r2;
        p = rand();
        if p < 0.5
            if abs(A) < 1
                % Encircling prey (exploitation around best candidate)
                D = abs(C * globalBest - population(i, :));
                newPopulation(i, :) = globalBest - A * D;
            else
                % Exploration via random search agent
                randIndex = randi([1 populationSize]);
                X_rand = population(randIndex, :);
                D = abs(C * X_rand - population(i, :));
                newPopulation(i, :) = X_rand - A * D;
            end
        else
            % Bubble-net attacking maneuver (spiral position update)
            D = abs(globalBest - population(i, :));
            newPopulation(i, :) = D * exp(-iter / maxIterations) * cos(2 * pi * iter / maxIterations) + globalBest;
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

save("WOA-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
