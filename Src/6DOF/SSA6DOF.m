function SSA6DOF(inputData)
% SSA6DOF  Salp Swarm Algorithm for 6-DOF aircraft parameter identification (60 derivatives).
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
foodPosition = population(1, :);
foodFitness = inf;

%% Optimization Loop
for iter = 1:maxIterations
    fitnessValues = EvaluateFitness6DOF(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    if localBestFitness < foodFitness
        foodPosition = population(idx, :);
        foodFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, :) = fitnessValues;
    globalrmse(iter) = mean(fitnessValues);

    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("SSA6DOF: convergence criterion reached");
        break;
    end

    % Leader-follower position update with exponential decay parameter c1
    c1 = 2 * exp(-(4 * iter / maxIterations)^2);
    for i = 1:populationSize
        if i <= populationSize / 2
            % Swarm leaders explore around current food source
            for j = 1:parameterCount
                c2 = rand(); c3 = rand();
                if c3 < 0.5
                    population(i, j) = foodPosition(j) + c1 * ((upperBounds(j) - lowerBounds(j)) * c2 + lowerBounds(j));
                else
                    population(i, j) = foodPosition(j) - c1 * ((upperBounds(j) - lowerBounds(j)) * c2 + lowerBounds(j));
                end
            end
        else
            % Followers follow preceding salp
            population(i, :) = 0.5 * (population(i, :) + population(i - 1, :));
        end
        population(i, :) = max(lowerBounds, min(upperBounds, population(i, :)));
    end
end

bestSolution = foodPosition;
elapsedTime = toc;
fprintf('SSA6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== SSA6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('SSA-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
