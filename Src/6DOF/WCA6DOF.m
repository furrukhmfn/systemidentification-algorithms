function WCA6DOF(inputData)
% WCA6DOF  Water Cycle Algorithm for 6-DOF aircraft parameter identification (60 derivatives).
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

dmax = 1e-5;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

population = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);

%% Optimization Loop
for iter = 1:maxIterations
    fitnessValues = EvaluateFitness6DOF(population, staticData);

    [sortedFitness, sortedIdx] = sort(fitnessValues);
    population = population(sortedIdx, :);
    sea = population(1, :);

    globalFitnessMatrix(iter, :) = fitnessValues;
    globalrmse(iter) = mean(fitnessValues);

    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("WCA6DOF: convergence criterion reached");
        break;
    end

    % Stream and river flow toward the best solution (sea) with evaporation-rain process
    C = 2;
    for i = 2:populationSize
        population(i, :) = population(i, :) + C * rand() * (sea - population(i, :));
        population(i, :) = max(lowerBounds, min(upperBounds, population(i, :)));

        if norm(population(i, :) - sea) < dmax
            population(i, :) = lowerBounds + rand(1, parameterCount) .* (upperBounds - lowerBounds);
        end
    end
end

bestSolution = sea;
elapsedTime = toc;
fprintf('WCA6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== WCA6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('WCA-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
