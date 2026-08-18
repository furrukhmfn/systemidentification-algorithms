function GWO6DOF(inputData)
% GWO6DOF  Grey Wolf Optimizer for full 6-DOF aircraft state-space identification.
% Estimates 60 parameters (18 stability + 42 control derivatives) across 12 rigid-body states:
%   States  (12): [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
%   Inputs   (7): [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]

arguments
    inputData InputData
end

%% Flight data and problem setup
staticData = PrepareFlightData6DOF(inputData);

% Unpack kinematics, Euler angles, and control inputs for post-processing
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
parameterCount = inputData.ParameterCount;   % 60 parameters
valueOfGravitationConstant = inputData.ValueOfGravitationConstant;

tic;
globalFitnessMatrix = zeros(maxIterations, populationSize);
globalrmse = [];

%% Initialize wolf population
wolves = InitializePopulation(populationSize, parameterCount, lowerBounds, upperBounds);

%% Main optimization loop
for iter = 1:maxIterations

    % Cost function evaluation across 6-DOF dynamics
    fitnessValues = EvaluateFitness6DOF(wolves, staticData);

    % Rank pack hierarchy: alpha (best), beta (2nd), delta (3rd)
    [fitnessValues, idx] = sort(fitnessValues);
    wolves = wolves(idx, :);

    alphaWolf = wolves(1, :);
    betaWolf  = wolves(2, :);
    deltaWolf = wolves(3, :);

    % Linearly decrease exploration factor
    a = 2 - iter * (2 / maxIterations);

    % Update agent positions based on leaders
    for i = 1:populationSize
        for j = 1:parameterCount

            % Alpha contribution
            r1 = rand();  r2 = rand();
            A1 = 2 * a * r1 - a;
            C1 = 2 * r2;
            D_alpha = abs(C1 * alphaWolf(j) - wolves(i, j));
            X1 = alphaWolf(j) - A1 * D_alpha;

            % Beta contribution
            r1 = rand();  r2 = rand();
            A2 = 2 * a * r1 - a;
            C2 = 2 * r2;
            D_beta = abs(C2 * betaWolf(j) - wolves(i, j));
            X2 = betaWolf(j) - A2 * D_beta;

            % Delta contribution
            r1 = rand();  r2 = rand();
            A3 = 2 * a * r1 - a;
            C3 = 2 * r2;
            D_delta = abs(C3 * deltaWolf(j) - wolves(i, j));
            X3 = deltaWolf(j) - A3 * D_delta;

            % Update position toward leadership centroid
            wolves(i, j) = (X1 + X2 + X3) / 3;

            % Enforce parameter search boundaries
            wolves(i, j) = max(lowerBounds(j), min(upperBounds(j), wolves(i, j)));
        end
    end

    globalFitnessMatrix(iter, :) = fitnessValues;
    globalrmse(iter) = mean(fitnessValues);

    % Relative convergence check on mean normalized fitness
    if iter > 1
        rel_diff = abs(globalrmse(iter) - globalrmse(iter - 1)) / (abs(globalrmse(iter - 1)) + eps);
        if rel_diff < 1e-5
            disp('GWO6DOF: convergence criterion reached');
            break;
        end
    end
end

globalBest   = alphaWolf;
bestSolution = globalBest;

elapsedTime = toc;
fprintf('GWO6DOF Execution time: %.5f seconds\n', elapsedTime);

%% Extract 6-DOF state-space matrices and export
u = staticData.u_all(end);
[A, B] = formatParameters6DOF(bestSolution, u);

disp('=== GWO6DOF: Best Solution Found ===');
fprintf('A Matrix (12x12):\n');  disp(A);
fprintf('B Matrix (12x7):\n');   disp(B);

save('GWO-6DOF-AlgoData.mat');
save(inputData.WorkspaceFileName);

end
