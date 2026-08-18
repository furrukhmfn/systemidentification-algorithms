function LSHADE6DOF(inputData)
% LSHADE6DOF  Linear Population Size Reduction SHADE for 6-DOF identification (60 derivatives).
% States (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
% Inputs  (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]

arguments
    inputData InputData
end

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

% History archive and adaptation memory configuration
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

populationSize = initialPopulationSize;
population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
globalBest = population(1, :);
globalBestFitness = inf;

%% Optimization Loop
for iter = 1:maxIterations
    % Linear population size reduction (LPSR)
    populationSize = round(initialPopulationSize - (iter / maxIterations) * (initialPopulationSize - minPopulationSize));
    population = population(1:populationSize, :);

    fitnessValues = EvaluateFitness6DOF(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    localBest = population(idx, :);
    if localBestFitness < globalBestFitness
        globalBest = localBest;
        globalBestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, 1:populationSize) = fitnessValues;

    globalrmse(iter) = mean(fitnessValues);
    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("LSHADE6DOF: convergence criterion reached");
        break;
    end

    memoryIndex = mod(iter - 1, memorySize) + 1;

    newPopulation = zeros(size(population));
    newFitness = fitnessValues;
    for i = 1:populationSize
        r1 = TournamentSelection(populationSize, 3, fitnessValues);
        r2 = TournamentSelection(populationSize, 3, fitnessValues);
        bestIdx = randi([1 populationSize]);

        % Adapt scaling factor and crossover rate from historical memory
        F(i) = MF(memoryIndex) + 0.1 * randn();
        CR(i) = MCR(memoryIndex) + 0.1 * randn();
        F(i) = min(max(F(i), 0.1), 0.9);
        CR(i) = min(max(CR(i), 0.1), 0.9);

        % Current-to-pbest mutation with external archive
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

        % Selection and archive update
        trialFitness = EvaluateFitness6DOF(trial, staticData);

        if trialFitness < newFitness(i)
            newPopulation(i, :) = trial;
            newFitness(i) = trialFitness;
            archive = [archive; population(i, :)];
        else
            newPopulation(i, :) = population(i, :);
        end
    end

    population = newPopulation;
end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('LSHADE6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== LSHADE6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('LSHADE-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end

function selected = TournamentSelection(popSize, tourSize, fitness)
    candidates = randi([1 popSize], 1, tourSize);
    [~, bestIdx] = min(fitness(candidates));
    selected = candidates(bestIdx);
end
