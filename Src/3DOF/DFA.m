function DFA(inputData)
% DFA  Dragonfly Algorithm for 24-parameter longitudinal aircraft model identification.
% Estimates aerodynamic stability and multi-control derivatives from flight data.

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

%% Initialize population and velocity vectors
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
deltaX = zeros(populationSize, parameterCount);
globalBest = population(1, :);
globalBestFitness = inf;

%% Main optimization loop
for iter = 1:maxIterations
    % Cost function evaluation
    fitnessValues = EvaluateFitness24(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx, :);
    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    % Track worst solution as repulsive enemy source
    [~, worstIdx] = max(fitnessValues);
    enemy = population(worstIdx, :);

    globalFitnessMatrix(iter, :) = fitnessValues;

    % Early convergence check
    globalrmse(iter) = mean(fitnessValues);
    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("algorithm end condition reached");
        break;
    end

    % Adaptive inertia weight and behavioral weighting factors
    w = 0.9 - iter * (0.9 - 0.4) / maxIterations;
    s = 0.1; a = 0.1; c = 0.1; f = 0.1; e = 0.1;

    % Swarm interaction: separation (S), alignment (A), and cohesion (C)
    newPopulation = zeros(populationSize, parameterCount);
    for i = 1:populationSize
        S = zeros(1, parameterCount);
        A = zeros(1, parameterCount);
        C = zeros(1, parameterCount);

        for j = 1:populationSize
            if i ~= j
                d = norm(population(i, :) - population(j, :));
                if d > 0
                    S = S - (population(i, :) - population(j, :)) / d;
                    A = A + population(j, :) / d;
                    C = C + population(j, :) / d;
                end
            end
        end

        % Food attraction (F) toward global best and enemy avoidance (E)
        F = globalBest - population(i, :);
        E = enemy + population(i, :);
        deltaX(i, :) = w * deltaX(i, :) + s * S + a * A + c * C + f * F + e * E;
        newPopulation(i, :) = population(i, :) + deltaX(i, :);
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

save("DFA-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
