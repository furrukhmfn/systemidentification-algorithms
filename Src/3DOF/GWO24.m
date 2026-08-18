function GWO24(inputData)
% GWO24  Grey Wolf Optimizer for 24-parameter longitudinal aircraft parameter estimation.
% Identifies stability and control derivatives using an alpha-beta-delta pack hierarchy.

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

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

%% Initialize wolf population
wolves = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);

%% Main optimization loop
for iter = 1:maxIterations
    % Cost function evaluation
    fitnessValues = EvaluateFitness24(wolves, staticData);

    % Rank pack members to establish leadership hierarchy
    [fitnessValues, idx] = sort(fitnessValues);
    wolves = wolves(idx, :);

    alphaWolf = wolves(1, :);
    betaWolf = wolves(2, :);
    deltaWolf = wolves(3, :);

    % Exploration-exploitation decay factor
    a = 2 - iter * (2 / maxIterations);
    for i = 1:populationSize
        for j = 1:parameterCount
            % Alpha wolf guidance
            r1 = rand();
            r2 = rand();
            A1 = 2 * a * r1 - a;
            C1 = 2 * r2;
            D_alpha = abs(C1 * alphaWolf(j) - wolves(i,j));
            X1 = alphaWolf(j) - A1 * D_alpha;

            % Beta wolf guidance
            r1 = rand();
            r2 = rand();
            A2 = 2 * a * r1 - a;
            C2 = 2 * r2;
            D_beta = abs(C2 * betaWolf(j) - wolves(i,j));
            X2 = betaWolf(j) - A2 * D_beta;

            % Delta wolf guidance
            r1 = rand();
            r2 = rand();
            A3 = 2 * a * r1 - a;
            C3 = 2 * r2;
            D_delta = abs(C3 * deltaWolf(j) - wolves(i,j));
            X3 = deltaWolf(j) - A3 * D_delta;

            % Update position toward leadership centroid
            wolves(i,j) = (X1 + X2 + X3) / 3;

            % Enforce parameter search boundaries
            wolves(i,j) = max(lowerBounds(j), min(upperBounds(j), wolves(i,j)));
        end
    end

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

globalBest = alphaWolf;
bestSolution = globalBest;
elapsedTime = toc;
fprintf('Execution time: %.5f seconds\n', elapsedTime);

%% Extract state-space model and save results
u = staticData.u_all(end);
[A, B] = formatParameters24(bestSolution, u);
disp('Best solution found. Matrices');
disp("A Matrix")
disp(A);
disp("B Matrix")
disp(B);

save("GWO-New-LongAlgoData.mat");
save(inputData.WorkspaceFileName);

end
