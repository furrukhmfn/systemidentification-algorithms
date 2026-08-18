function GHOA6DOF(inputData)
% GHOA6DOF  Grasshopper Optimization Algorithm for 6-DOF aircraft identification (60 derivatives).
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

cMin = 0.004;
cMax = 1;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);

globalBest = population(1, :);
globalBestFitness = inf;

%% Optimization Loop
for iter = 1:maxIterations
    % Decreasing comfort-zone coefficient balancing exploration and exploitation
    c = cMax - iter * (cMax - cMin) / maxIterations;

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
        disp("GHOA6DOF: convergence criterion reached");
        break;
    end

    % Social interaction and attraction toward target solution
    for i = 1:populationSize
        S_i = zeros(1, parameterCount);
        for j = 1:populationSize
            if i ~= j
                dist = norm(population(j, :) - population(i, :));
                r_ij_vec = (population(j, :) - population(i, :)) / (dist + eps);
                s_ij = 0.5 * exp(-dist / 1.5) - exp(-dist);
                S_i = S_i + c * ((upperBounds - lowerBounds) / 2) * s_ij * r_ij_vec;
            end
        end
        X_new = c * S_i + globalBest;
        population(i, :) = max(lowerBounds, min(upperBounds, X_new));
    end
end

bestSolution = globalBest;
elapsedTime = toc;
fprintf('GHOA6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== GHOA6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('GHOA-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
