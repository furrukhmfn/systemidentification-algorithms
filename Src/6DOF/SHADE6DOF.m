function SHADE6DOF(inputData)
% SHADE6DOF  Success-History Adaptive Differential Evolution for 6-DOF identification (60 derivatives).
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

% Historical parameter memories
H = 100;
MF = 0.5 * ones(H, 1);
MCR = 0.5 * ones(H, 1);
k = 1;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

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
    globalrmse(iter) = mean(fitnessValues);

    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("SHADE6DOF: convergence criterion reached");
        break;
    end

    newPopulation = zeros(populationSize, parameterCount);
    for i = 1:populationSize
        idx_m = randi([1 H]);
        CR = normrnd(MCR(idx_m), 0.1);
        CR = min(max(CR, 0), 1);

        % Cauchy-distributed mutation scale factor
        F = cauchyrnd(MF(idx_m), 0.1);
        while F <= 0
            F = cauchyrnd(MF(idx_m), 0.1);
        end
        F = min(F, 1);

        r1 = randi([1 populationSize]);
        r2 = randi([1 populationSize]);

        % Current-to-best mutation
        mutant = population(i, :) + F * (globalBest - population(i, :)) + F * (population(r1, :) - population(r2, :));
        mutant = max(lowerBounds, min(upperBounds, mutant));

        % Binomial crossover
        trial = population(i, :);
        j_rand = randi([1 parameterCount]);
        for j = 1:parameterCount
            if rand() <= CR || j == j_rand
                trial(j) = mutant(j);
            end
        end
        newPopulation(i, :) = trial;
    end
    population = newPopulation;
end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('SHADE6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== SHADE6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('SHADE-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end

function val = cauchyrnd(location, scale)
    val = location + scale * tan(pi * (rand() - 0.5));
end
