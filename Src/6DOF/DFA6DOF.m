function DFA6DOF(inputData)
% DFA6DOF  Dragonfly Algorithm for full 6-DOF aircraft parameter identification (60 derivatives).
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

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
deltaX = zeros(populationSize, parameterCount);
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

    [~, worstIdx] = max(fitnessValues);
    enemy = population(worstIdx, :);

    globalFitnessMatrix(iter, :) = fitnessValues;
    globalrmse(iter) = mean(fitnessValues);

    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("DFA6DOF: convergence criterion reached");
        break;
    end

    % Dynamic swarming weights (separation, alignment, cohesion, food attraction, enemy distraction)
    w = 0.9 - iter * (0.9 - 0.4) / maxIterations;
    s = 0.1; a = 0.1; c = 0.1; f = 0.1; e = 0.1;

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
fprintf('DFA6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== DFA6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('DFA-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
