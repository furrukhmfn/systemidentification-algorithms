function WOA6DOF(inputData)
% WOA6DOF  Whale Optimization Algorithm for 6-DOF aircraft parameter identification (60 derivatives).
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
bestPosition = population(1, :);
bestFitness = inf;

%% Optimization Loop
for iter = 1:maxIterations
    fitnessValues = EvaluateFitness6DOF(population, staticData);

    [localBestFitness, idx] = min(fitnessValues);
    if localBestFitness < bestFitness
        bestPosition = population(idx, :);
        bestFitness = localBestFitness;
    end

    globalFitnessMatrix(iter, :) = fitnessValues;
    globalrmse(iter) = mean(fitnessValues);

    if iter > 1 && abs(globalrmse(iter) - globalrmse(iter - 1)) < 0.0001
        disp("WOA6DOF: convergence criterion reached");
        break;
    end

    % Encircling prey and bubble-net spiral hunting mechanisms
    a = 2 - iter * (2 / maxIterations);
    for i = 1:populationSize
        r1 = rand(); r2 = rand();
        A_val = 2 * a * r1 - a;
        C_val = 2 * r2;
        p_prob = rand();
        l = (rand() - 0.5) * 2;

        if p_prob < 0.5
            if abs(A_val) < 1
                % Encircling best solution
                D = abs(C_val * bestPosition - population(i, :));
                population(i, :) = bestPosition - A_val * D;
            else
                % Exploration via random search agent
                rand_idx = randi([1 populationSize]);
                X_rand = population(rand_idx, :);
                D = abs(C_val * X_rand - population(i, :));
                population(i, :) = X_rand - A_val * D;
            end
        else
            % Spiral bubble-net trajectory
            D1 = abs(bestPosition - population(i, :));
            population(i, :) = D1 * exp(b * l) * cos(2 * pi * l) + bestPosition;
        end
        population(i, :) = max(lowerBounds, min(upperBounds, population(i, :)));
    end
end

bestSolution = bestPosition;
elapsedTime = toc;
fprintf('WOA6DOF Execution time: %.5f seconds\n', elapsedTime);

% Format continuous-time state-space matrices
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== WOA6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('WOA-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
