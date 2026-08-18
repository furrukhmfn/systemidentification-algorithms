function GHOA24(inputData)
% GHOA24  Grasshopper Optimization Algorithm for 24-parameter longitudinal model identification.
% Estimates aircraft stability and multi-surface control derivatives from flight test data.

arguments
    inputData InputData
end

%% Flight data and problem setup
staticData = PrepareFlightData24(inputData);

% Unpack raw data into local workspace variables for post-processing
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

% Comfort zone contraction limits
cMin = 0.004;
cMax = 1;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

%% Initialize population and best solution tracking
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
globalBest = population(1,:);
globalBestFitness = inf;

%% Main optimization loop
for iter = 1:maxIterations
    % Linearly decrease comfort zone parameter
    c = cMax - iter * (cMax - cMin) / maxIterations;
    
    % Evaluate candidate model fitness
    fitnessValues = EvaluateFitness24(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx, :);

    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    % Update agent positions via social forces, gravity, and wind drift
    newPopulation = zeros(populationSize, parameterCount);

    for i = 1:populationSize
        S = zeros(1, parameterCount);  % Social interaction vector

        % Attraction/repulsion forces between agents
        for j = 1:populationSize
            if i ~= j
                distance = norm(population(j,:) - population(i,:));

                direction = (population(j,:) - population(i,:)) / ...
                            (distance + eps);

                S = S + (c / 2) .* (upperBounds - lowerBounds) .* ...
                        direction .* exp(distance / 1.5);
            end
        end

        gravity = -c * globalBest;                   % Attraction toward global best
        wind = c * rand(1, parameterCount);         % Stochastic drift component

        newPosition = population(i,:) + S + gravity + wind;

        % Enforce parameter search boundaries
        newPopulation(i, :) = max(lowerBounds, ...
                                min(upperBounds, newPosition));
    end

    population = newPopulation;
    globalFitnessMatrix(iter, :) = fitnessValues;

    % Early convergence check
    globalrmse(iter) = mean(fitnessValues);
    if iter > 1
        diff_val = abs(globalrmse(iter) - globalrmse(iter - 1));
        if diff_val < 0.0001
            disp("algorithm end condition reached");
            break;
        end
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

save("GHOA-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
