function LSHADE(inputData)
% LSHADE  Linear Success-History Based Adaptive Differential Evolution for 24-parameter longitudinal model identification.
% Features linear population size reduction (LPSR) and historical memory adaptation of mutation and crossover rates.

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

% LSHADE hyperparameter memory and population limits
initialPopulationSize = populationSize;
minPopulationSize = 100;
memorySize = 100;
F = rand(initialPopulationSize, 1);
CR = rand(initialPopulationSize, 1);
archive = [];
MF = rand(memorySize, 1);
MCR = rand(memorySize, 1);

tic;
globalFitnessMatrix = zeros(maxIterations, initialPopulationSize);
globalrmse = [];

%% Initialize population
populationSize = initialPopulationSize;
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
globalBest = population(1, :);
globalBestFitness = inf;

%% Main optimization loop
for iter = 1:maxIterations
    % Linear population size reduction (LPSR)
    populationSize = round(initialPopulationSize - (iter / maxIterations) * (initialPopulationSize - minPopulationSize));
    population = population(1:populationSize, :);

    % Fitness evaluation
    fitnessValues = EvaluateFitness24(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx, :);
    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, 1:populationSize) = fitnessValues;

    % Early convergence check
    globalrmse(iter) = mean(fitnessValues);
    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("algorithm end condition reached");
        break;
    end

    memoryIndex = mod(iter - 1, memorySize) + 1;

    % Mutation and crossover with historical parameter sampling
    newPopulation = zeros(size(population));
    newFitness = fitnessValues;
    for i = 1:populationSize
        r1 = TournamentSelection(populationSize, 3, fitnessValues);
        r2 = TournamentSelection(populationSize, 3, fitnessValues);
        bestIdx = randi([1 populationSize]);

        % Sample F and CR from historical memory
        F(i) = MF(memoryIndex) + 0.1 * randn();
        CR(i) = MCR(memoryIndex) + 0.1 * randn();
        F(i) = min(max(F(i), 0.1), 0.9);
        CR(i) = min(max(CR(i), 0.1), 0.9);

        % Current-to-pbest mutation strategy
        mutant = population(i, :) + F(i) * (population(bestIdx, :) - population(i, :)) + ...
                 F(i) * (population(r1, :) - population(r2, :));
        mutant = max(lowerBounds, min(upperBounds, mutant));

        % Binomial crossover
        trial = population(i, :);
        for j = 1:parameterCount
            if rand() <= CR(i) || j == randi([1 parameterCount])
                trial(j) = mutant(j);
            end
        end
        trial = max(lowerBounds, min(upperBounds, trial));

        trialFitness = EvaluateFitness24(trial, staticData);

        % Selection and external archive update
        if trialFitness < newFitness(i)
            newPopulation(i, :) = trial;
            newFitness(i) = trialFitness;
            archive = [archive; population(i, :)];
        else
            newPopulation(i, :) = population(i, :);
        end
    end

    population = newPopulation;
    fitnessValues = newFitness;

    % Adapt historical memory means
    if sum(fitnessValues < fitnessValues) > 0
        MF(memoryIndex) = mean(F(fitnessValues < fitnessValues));
        MCR(memoryIndex) = mean(CR(fitnessValues < fitnessValues));
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

save("LSHADE-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end

function idx = TournamentSelection(popSize, k, fitness)
    candidates = randi(popSize, k, 1);
    [~, best] = min(fitness(candidates));
    idx = candidates(best);
end
