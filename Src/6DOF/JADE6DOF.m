function JADE6DOF(inputData)
% JADE6DOF  Adaptive Differential Evolution for 6-DOF aircraft parameter identification (60 derivatives).
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

b = 1;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);
bestFlames = population;

%% Optimization Loop
for iter = 1:maxIterations
    fitnessValues = EvaluateFitness6DOF(population, staticData);

    [sortedFitnessValues, sortedIndices] = sort(fitnessValues);
    bestFlames = population(sortedIndices, :);

    globalFitnessMatrix(iter, :) = fitnessValues;
    globalrmse(iter) = mean(fitnessValues);

    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("JADE6DOF: convergence criterion reached");
        break;
    end

    % Spiral search update around sorted elite solutions
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
fprintf('JADE6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== JADE6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('JADE-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
