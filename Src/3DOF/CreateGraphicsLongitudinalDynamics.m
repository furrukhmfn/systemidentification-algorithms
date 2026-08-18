function CreateGraphicsLongitudinalDynamics(inputGraphicsData)
% CreateGraphicsLongitudinalDynamics Generates validation plots (RMSE progression,
% derivative tracking, step response, and convergence curves) for identified
% 15-parameter longitudinal dynamics models.
arguments
    inputGraphicsData InputGraphicsData;
end

%% Load Identified Workspace
load(inputGraphicsData.OutputDataFileName);

algorithmName = inputGraphicsData.AlgorithmName;
caseNumber = inputGraphicsData.CaseNumber;

%% RMSE Convergence Plots
fig = figure;
i = 1:populationSize;

subplot(2,2,1)
j=[1,2,3,4,5];
j=j(j <= size(globalFitnessMatrix,1));
if ~isempty(j)
    plot(i,globalFitnessMatrix(j,:)');
    legend(cellstr(num2str(j', 'it-%d')), 'FontName', 'Times New Roman');
end

subplot(2,2,2)
j=[6,7,8,9,10];
j=j(j <= size(globalFitnessMatrix,1));
if ~isempty(j)
    plot(i,globalFitnessMatrix(j,:)');
    legend(cellstr(num2str(j', 'it-%d')), 'FontName', 'Times New Roman');
end

subplot(2,2,3)
j=[11,12,13,14,15];
j=j(j <= size(globalFitnessMatrix,1));
if ~isempty(j)
    plot(i,globalFitnessMatrix(j,:)');
    legend(cellstr(num2str(j', 'it-%d')), 'FontName', 'Times New Roman');
end

subplot(2,2,4)
j=[16,17,18,19,20];
j=j(j <= size(globalFitnessMatrix,1));
if ~isempty(j)
    plot(i,globalFitnessMatrix(j,:)');
    legend(cellstr(num2str(j', 'it-%d')), 'FontName', 'Times New Roman');
end

han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'RMSE', 'FontName', 'Times New Roman');
xlabel(han,'Iterations', 'FontName', 'Times New Roman');
mainTitle = sgtitle(algorithmName + " - " + caseNaming(caseNumber) + " - Longitudinal", 'FontName', 'Times New Roman');
mainTitle.FontSize = 10;
fontname(fig, 'Times New Roman');
saveas(gcf, "output/Error" + algorithmName + "_" + "Case" + num2str(caseNumber) + "_Longitudinal.png");

%% Derivative State Tracking Validation
thrust_val = 0;

u_all = Vb(:, 1);
w_all = Vb(:, 3);
q_all = pqr(:, 2);
theta_all = phi_theta_psi(2, :)';
h_all = -Xe(:, 3);
elevator_all = elevator(:, 1);

u_dot_all = Accels(:, 1);
w_dot_all = Accels(:, 3);
q_dot_all = pdot_qdot_rdot(:, 2);
h_dot_all = Ve(:, 3);
phi_all = phi_theta_psi(1, :)';

% Kinematic Euler pitch rate
theta_dot_all = cos(phi_all) .* q_all - sin(phi_all) .* pqr(:, 3);

p1 = globalBest(1);
p2 = globalBest(2);
p3 = globalBest(3);
p4 = globalBest(4);
p5 = globalBest(5);
p6 = globalBest(6);
p7 = globalBest(7);
p8 = globalBest(8);
p9 = globalBest(9);
p10 = globalBest(10);
p11 = globalBest(11);
p12 = globalBest(12);
p13 = globalBest(13);
p14 = globalBest(14);
p15 = globalBest(15);

xdotcheck_1 = p1 .* u_all + p2 .* w_all + p3 .* q_all - valueOfGravitationConstant .* theta_all + p10 .* elevator_all + p11 .* thrust_val;
xdotcheck_2 = p4 .* u_all + p5 .* w_all + p6 .* q_all + p12 .* elevator_all + p13 .* thrust_val;
xdotcheck_3 = p7 .* u_all + p8 .* w_all + p9 .* q_all + p14 .* elevator_all + p15 .* thrust_val;
xdotcheck_4 = q_all;
xdotcheck_5 = -w_all + u_all .* theta_all;

xdotcheck = [xdotcheck_1, xdotcheck_2, xdotcheck_3, xdotcheck_4, xdotcheck_5]';
xdotreal = [u_dot_all, w_dot_all, q_dot_all, theta_dot_all, h_dot_all]';
it = size(Vb, 1);

u = u_all(end);
[A, B] = formatParameters(globalBest, u, valueOfGravitationConstant);

fig4 = figure;
subplot(3,2,1);
plot(1:it, xdotreal(1,:)); hold on;
plot(1:it, xdotcheck(1,:));
title('$\dot{u}$', 'Interpreter','latex');
legend({'real', 'sim'}, 'FontName', 'Times New Roman');
ylabel('m/s^2');

subplot(3,2,2);
plot(1:it, xdotreal(2,:)); hold on;
plot(1:it, xdotcheck(2,:));
legend({'real', 'sim'});
title('$\dot{w}$', 'Interpreter','latex');
ylabel('m/s^2');

subplot(3,2,3);
plot(1:it, xdotreal(3,:)); hold on;
plot(1:it, xdotcheck(3,:));
legend({'real', 'sim'});
title('$\dot{q}$', 'Interpreter','latex');
ylabel('rad/s^2');

subplot(3,2,4);
plot(1:it, xdotreal(4,:)); hold on;
plot(1:it, xdotcheck(4,:));
legend({'real', 'sim'});
title('$\dot{\theta}$', 'Interpreter','latex');
ylabel('rad/s');

subplot(3,2,[5,6]);
plot(1:it, xdotreal(5,:)); hold on;
plot(1:it, -xdotcheck(5,:));
xlabel('Iteration');
legend({'real', 'sim'});
title('$\dot{h}$', 'Interpreter','latex');
ylabel('m/s');
mainTitle = sgtitle("Validation - "+algorithmName + " - " + caseNaming(caseNumber) + " - Longitudinal");
mainTitle.FontSize = 10;
fontname(fig4, 'Times New Roman');
saveas(gcf, "Validation_" + algorithmName + "_" + "Case" + num2str(caseNumber) + "_Longitudinal.png");

%% Step Response Analysis
C = [eye(5)];
C(5,:) = [0 1/u 0 0 0];
C = -C;
D = -ones(2,5)';
sys= ss(A,B,C,D);
P  = tf(sys)
predicted = step(P,1);

% Reference aircraft model
AReal = [-0.0013         0.012   1.1229    9.8100         0
   -0.2363    0.9317    0.3597         0         0
         0    0.6137    1.4657         0         0
         0         0    1.0000         0         0
         0   -1.0000         0  114.6612         0];

BReal = [ 0.0472    0.7602
    0.0864    0.5424
   -0.6128         0
         0         0
         0         0];

sysReal= ss(AReal,BReal,C,D);
PReal  = tf(sysReal)
predictedReal = step(PReal,1);

figure;
h = gobjects(2, 1); 

for varIdx = 1:5
    for inputIdx = 1:2
        subplotIdx = (varIdx - 1) * 2 + inputIdx;
        subplot(5, 2, subplotIdx);

        realData = squeeze(predictedReal(:, varIdx, inputIdx));
        predictedData = squeeze(predicted(:, varIdx, inputIdx));

        h(1) = stairs(realData, 'b', 'LineWidth', 1.5); hold on;
        h(2) = stairs(predictedData, 'r--', 'LineWidth', 1.5); hold off;

       set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
        grid on;
    end
end

han = axes(gcf, 'visible', 'off'); 
han.XLabel.Visible = 'on'; 
xlabel(han, 'Time Step', 'FontName', 'Times New Roman', 'FontSize', 12);

sgtitle("Step Plot" + " - " + caseNaming(caseNumber), 'FontName', 'Times New Roman', 'FontSize', 12);

legend(h, {'Real', 'Predicted'}, 'Position', [0.8 0.90 0.1 0.05], ...
    'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, "Step" + algorithmName + "_" + "Case" + num2str(caseNumber) + "_Longitudinal.png");

%% Learning Curve and Statistical Export
meanFitness = mean(globalFitnessMatrix');

fig = figure('Color', 'w'); 
hold on; grid on; box on;

plot(1: iter, meanFitness, '-k', 'LineWidth', 1.5, 'DisplayName', 'sin(x)');

xlabel('Iteration', 'FontName', 'Times New Roman');
ylabel('Mean Error', 'FontName', 'Times New Roman');
title("Mean Error - "+ algorithmName + " - " + caseNaming(caseNumber) + " - Longitudinal.png", 'FontName', 'Times New Roman');

ax = gca;
ax.Position = [0.12 0.12 0.8 0.8];
fontname(fig, 'Times New Roman');
f = gcf;
exportgraphics(f,"MeanErrorGraph"+ algorithmName + "_" + "Case" + num2str(caseNumber) + "_Longitudinal.png",'Resolution',600)

% Export parameters and Friedman test ranking to Excel
excelSheetFileName = "Data_"+"Case" + num2str(caseNumber) + "_Longitudinal.xlsx";

try
    globalFitnessMatrix = globalFitnessMatrix'; 
    [p, tbl, stats] = friedman(globalFitnessMatrix, iter, 'off');
    
    xlswrite(excelSheetFileName, bestSolution, 'Sheet1');
    tblArray = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
    writetable(tblArray, excelSheetFileName, 'Sheet', 'Sheet2');
catch ME
    warning('Could not perform Friedman test or write Sheet2. Error: %s', ME.message);
    try
        xlswrite(excelSheetFileName, bestSolution, 'Sheet1');
    catch xlME
        warning('Failed to write Sheet1: %s', xlME.message);
    end
end

end
