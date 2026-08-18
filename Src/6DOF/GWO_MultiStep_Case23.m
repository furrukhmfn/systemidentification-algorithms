function GWO_MultiStep_Case23()
% GWO_MultiStep_Case23  Comparative analysis of Single-Step vs Multi-Step GWO (Case 23).
% Workflow:
%   Pass 1 (Single-Step): Global search across baseline physical bounds.
%   Pass 2 (Multi-Step) : Refined local search within narrowed envelope around Pass 1 solution.
% Tracks per-iteration execution timing, ODE trajectory matching, and statistical metrics.

%% Configuration
CASE_NUM=23; POP_SIZE=6000; OFFSET_INC=0; MAX_ITER=20; PARAM_COUNT=24;
SHRINK_FACTOR=0.30; VALIDATION_MARGIN=1000;
accelFile='cruise_acceleration.mat'; statesFile='cruise_outputStates.mat'; accutFile='cruise_acctuators.mat';
if ~exist('output','dir'), mkdir('output'); end
outBase=fullfile('output',sprintf('GWO_MultiStep_Case%d',CASE_NUM));
diary([outBase '_log.txt']);
algoDir=fullfile(fileparts(mfilename('fullpath')),'more_algorithms');
if exist(algoDir,'dir'), addpath(algoDir); end
lb_orig=[-10,-100,-10,-10,-10,-10,-1,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20,-20];
ub_orig=[0,100,10,0,-eps,10,1,0,0,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20];
workFile=fullfile('output',sprintf('GWO_Case%d_24.mat',CASE_NUM));
inputData=InputData(accelFile,statesFile,accutFile,workFile,lb_orig,ub_orig,MAX_ITER,POP_SIZE,PARAM_COUNT);
if OFFSET_INC~=0, inputData=inputData.IncreamentOffSet(OFFSET_INC); end
fprintf('\n================================================================\n');
fprintf('  GWO Multi-Step vs Single-Step -- Case (%d,%d)  PopSize=%d  MaxIter=%d  Shrink=%.0f%%\n',...
    floor(CASE_NUM/10),mod(CASE_NUM,10),POP_SIZE,MAX_ITER,SHRINK_FACTOR*100);
fprintf('================================================================\n\n');
staticData=PrepareFlightData24(inputData);

%% Pass 1: Global Search with Baseline Bounds
fprintf('\n--- PASS 1: Single-Step GWO (Original Bounds) ---\n');
[globalBest_p1,fitnessHistory_p1,convIter_p1,timePerIter_p1]=runGWO(staticData,lb_orig,ub_orig,POP_SIZE,MAX_ITER,PARAM_COUNT);
elapsedP1=sum(timePerIter_p1);
fit_p1=EvaluateFitness24(globalBest_p1,staticData);
fprintf('  Single-Step: %.2fs total | fitness=%.6f | conv_iter=%d/%d\n',elapsedP1,fit_p1,convIter_p1,numel(fitnessHistory_p1));
fprintf('  Per-iter: mean=%.2fs  min=%.2fs  max=%.2fs\n',mean(timePerIter_p1),min(timePerIter_p1),max(timePerIter_p1));

%% Refine Parameter Search Envelope
range_orig=ub_orig-lb_orig;
lb_p2=max(lb_orig,globalBest_p1-SHRINK_FACTOR*range_orig);
ub_p2=min(ub_orig,globalBest_p1+SHRINK_FACTOR*range_orig);
bad=lb_p2>=ub_p2; lb_p2(bad)=lb_orig(bad); ub_p2(bad)=ub_orig(bad);
fprintf('\n  Bound width: %.4f -> %.4f\n',mean(ub_orig-lb_orig),mean(ub_p2-lb_p2));

%% Pass 2: Refined Local Search
fprintf('\n--- PASS 2: Multi-Step GWO (Refined Bounds) ---\n');
[globalBest_p2,fitnessHistory_p2,convIter_p2,timePerIter_p2]=runGWO(staticData,lb_p2,ub_p2,POP_SIZE,MAX_ITER,PARAM_COUNT);
elapsedP2=sum(timePerIter_p2);
fit_p2=EvaluateFitness24(globalBest_p2,staticData);
fprintf('  Multi-Step: %.2fs total | fitness=%.6f | conv_iter=%d/%d\n',elapsedP2,fit_p2,convIter_p2,numel(fitnessHistory_p2));
fprintf('  Per-iter: mean=%.2fs  min=%.2fs  max=%.2fs\n',mean(timePerIter_p2),min(timePerIter_p2),max(timePerIter_p2));
fprintf('  Grand total: %.2fs (%.2f min)\n',elapsedP1+elapsedP2,(elapsedP1+elapsedP2)/60);

%% ODE Trajectory Validation (ODE15s)
fprintf('\n--- ODE Validation ---\n');
Vb=staticData.Vb; pqr_mat=staticData.pqr; phi_theta_psi=staticData.phi_theta_psi;
Xe=staticData.Xe; elevator=staticData.elevator;
flapPos=staticData.flapPos; flapNeg=staticData.flapNeg; flapDiff=staticData.flapDiff;
time=staticData.time_all;
N=size(Vb,1); dS=1+VALIDATION_MARGIN; dE=N-VALIDATION_MARGIN;
u_data=Vb(dS:dE,1); w_data=Vb(dS:dE,3); q_data=pqr_mat(dS:dE,2);
theta_data=phi_theta_psi(2,dS:dE)'; h_data=-Xe(dS:dE,3);
de_data=elevator(dS:dE,1); fp_data=flapPos(dS:dE,1);
fn_data=flapNeg(dS:dE,1); fd_data=flapDiff(dS:dE,1);
X_real=[u_data,w_data,q_data,theta_data,h_data]'; X0=X_real(:,1);
thrust=0; uo_mean=mean(u_data);
fprintf('  Simulating Single-Step GWO (Pass 1)...\n');
X_sim_p1=simulateODE(globalBest_p1,uo_mean,time,u_data,de_data,fp_data,fn_data,fd_data,thrust,X0);
fprintf('  Simulating Multi-Step GWO (Pass 2)...\n');
X_sim_p2=simulateODE(globalBest_p2,uo_mean,time,u_data,de_data,fp_data,fn_data,fd_data,thrust,X0);

%% Error and Fit Metrics
stateNames={'u (m/s)','w (m/s)','q (rad/s)','theta (rad)','h (m)'}; nStates=5;
metricsP1=zeros(nStates,3); metricsP2=zeros(nStates,3);
for s=1:nStates
    y=X_real(s,:)';
    metricsP1(s,:)=computeMetrics3(y,X_sim_p1(:,s));
    metricsP2(s,:)=computeMetrics3(y,X_sim_p2(:,s));
end
convFrac_p1=convIter_p1/numel(fitnessHistory_p1);
convFrac_p2=convIter_p2/numel(fitnessHistory_p2);
timeToConv_p1=sum(timePerIter_p1(1:convIter_p1));
timeToConv_p2=sum(timePerIter_p2(1:convIter_p2));

agg.P1.MAE=mean(metricsP1(:,2)); agg.P1.RMSE=mean(metricsP1(:,1)); agg.P1.R2=mean(metricsP1(:,3));
agg.P1.Stability=std(metricsP1(:,1)); agg.P1.ConvSpeed=convFrac_p1;
agg.P1.TotalTime_s=elapsedP1; agg.P1.MeanIterTime=mean(timePerIter_p1);
agg.P1.MinIterTime=min(timePerIter_p1); agg.P1.MaxIterTime=max(timePerIter_p1);
agg.P1.TimeToConv_s=timeToConv_p1;

agg.P2.MAE=mean(metricsP2(:,2)); agg.P2.RMSE=mean(metricsP2(:,1)); agg.P2.R2=mean(metricsP2(:,3));
agg.P2.Stability=std(metricsP2(:,1)); agg.P2.ConvSpeed=convFrac_p2;
agg.P2.TotalTime_s=elapsedP2; agg.P2.MeanIterTime=mean(timePerIter_p2);
agg.P2.MinIterTime=min(timePerIter_p2); agg.P2.MaxIterTime=max(timePerIter_p2);
agg.P2.TimeToConv_s=timeToConv_p2;

%% Console Summary Report
caseDisplay=sprintf('(%d,%d)',floor(CASE_NUM/10),mod(CASE_NUM,10));
fprintf('\n================================================================\n');
fprintf('  METRICS SUMMARY -- Case %s\n',caseDisplay);
fprintf('  %-24s  %-9s  %-9s  %-9s  %-9s  %-7s\n','Model','MAE','RMSE','Stability','ConvFrac','R2');
fprintf('  %s\n',repmat('-',1,70));
fprintf('  %-24s  %-9.4f  %-9.4f  %-9.4f  %-9.4f  %-7.4f\n','Single-Step GWO (Pass-1)',agg.P1.MAE,agg.P1.RMSE,agg.P1.Stability,convFrac_p1,agg.P1.R2);
fprintf('  %-24s  %-9.4f  %-9.4f  %-9.4f  %-9.4f  %-7.4f\n','Multi-Step GWO (Pass-2)',agg.P2.MAE,agg.P2.RMSE,agg.P2.Stability,convFrac_p2,agg.P2.R2);
fprintf('================================================================\n\n');
fprintf('  TIMING SUMMARY:\n');
fprintf('  %-24s  %-11s  %-13s  %-12s  %-12s  %-13s\n','Model','Total(s)','MeanIter(s)','MinIter(s)','MaxIter(s)','TimeToConv(s)');
fprintf('  %s\n',repmat('-',1,90));
fprintf('  %-24s  %-11.2f  %-13.2f  %-12.2f  %-12.2f  %-13.2f\n','Single-Step GWO (Pass-1)',agg.P1.TotalTime_s,agg.P1.MeanIterTime,agg.P1.MinIterTime,agg.P1.MaxIterTime,agg.P1.TimeToConv_s);
fprintf('  %-24s  %-11.2f  %-13.2f  %-12.2f  %-12.2f  %-13.2f\n','Multi-Step GWO (Pass-2)',agg.P2.TotalTime_s,agg.P2.MeanIterTime,agg.P2.MinIterTime,agg.P2.MaxIterTime,agg.P2.TimeToConv_s);
fprintf('  GRAND TOTAL (Pass1+Pass2): %.2fs (%.2f min)\n',elapsedP1+elapsedP2,(elapsedP1+elapsedP2)/60);
fprintf('================================================================\n\n');

fprintf('  Per-State Detail:\n');
fprintf('  %-13s  %-12s  %-12s  %-10s  %-10s\n','State','Single RMSE','Multi RMSE','Single R2','Multi R2');
fprintf('  %s\n',repmat('-',1,62));
for s=1:nStates
    fprintf('  %-13s  %-12.4f  %-12.4f  %-10.4f  %-10.4f\n',...
        stateNames{s},metricsP1(s,1),metricsP2(s,1),metricsP1(s,3),metricsP2(s,3));
end

%% Validation Figures
applyIEEEDefaults();
cFlight=[0.10 0.10 0.10]; cPass1=[0.00 0.45 0.70]; cPass2=[0.80 0.20 0.10];
modelLabels={'Single-Step GWO','Multi-Step GWO'}; barColors=[cPass1;cPass2];

% Figure 1: State Validation Comparisons
fig1=figure('Units','inches','Position',[0.5 0.5 7.16 9.0],'Color','w');
for s=1:nStates
    subplot(3,2,s);
    plot(time,X_real(s,:),'-','Color',cFlight,'LineWidth',0.9,'DisplayName','Flight Data'); hold on;
    plot(time,X_sim_p1(:,s),'--','Color',cPass1,'LineWidth',0.8,'DisplayName','Single-Step GWO');
    plot(time,X_sim_p2(:,s),'-.','Color',cPass2,'LineWidth',0.8,'DisplayName','Multi-Step GWO');
    xlabel('Time (s)','FontSize',8); ylabel(stateNames{s},'FontSize',8);
    if s==1, legend('Location','best','FontSize',7); end
    set(gca,'FontName','Times New Roman','FontSize',8); grid on;
end
sgtitle(sprintf('State Validation -- Case %s',caseDisplay),'FontWeight','normal','FontSize',10,'FontName','Times New Roman');
saveFig(fig1,[outBase '_States.png'],[outBase '_States.eps']);
fprintf('  Saved: _States.png\n');

% Figure 2: Simulation Error Residuals
fig2=figure('Units','inches','Position',[0.5 0.5 7.16 9.0],'Color','w');
for s=1:nStates
    subplot(3,2,s);
    plot(time,X_real(s,:)'-X_sim_p1(:,s),'-','Color',cPass1,'LineWidth',0.7,'DisplayName','Single-Step Err'); hold on;
    plot(time,X_real(s,:)'-X_sim_p2(:,s),'-','Color',cPass2,'LineWidth',0.7,'DisplayName','Multi-Step Err');
    yline(0,'k--','LineWidth',0.5);
    xlabel('Time (s)','FontSize',8); ylabel(['Err: ' stateNames{s}],'FontSize',8);
    if s==1, legend('Location','best','FontSize',7); end
    set(gca,'FontName','Times New Roman','FontSize',8); grid on;
end
sgtitle(sprintf('Validation Error -- Case %s',caseDisplay),'FontWeight','normal','FontSize',10,'FontName','Times New Roman');
saveFig(fig2,[outBase '_Error.png'],[outBase '_Error.eps']);
fprintf('  Saved: _Error.png\n');

% Figure 3: Convergence History & Per-Iteration Computational Time
fig3=figure('Units','inches','Position',[0.5 0.5 7.16 3.0],'Color','w');
iP1=1:numel(fitnessHistory_p1); iP2=1:numel(fitnessHistory_p2);
subplot(1,2,1);
plot(iP1,fitnessHistory_p1,'-','Color',cPass1,'LineWidth',1.1,'DisplayName','Single-Step'); hold on;
plot(iP2,fitnessHistory_p2,'--','Color',cPass2,'LineWidth',1.1,'DisplayName','Multi-Step');
if convIter_p1<=numel(fitnessHistory_p1)
    plot(convIter_p1,fitnessHistory_p1(convIter_p1),'o','Color',cPass1,'MarkerSize',6,'LineWidth',1.2,'HandleVisibility','off');
    text(convIter_p1,fitnessHistory_p1(convIter_p1),sprintf('  i%d\n  %.1fs',convIter_p1,timeToConv_p1),'FontSize',6,'Color',cPass1,'FontName','Times New Roman');
end
if convIter_p2<=numel(fitnessHistory_p2)
    plot(convIter_p2,fitnessHistory_p2(convIter_p2),'s','Color',cPass2,'MarkerSize',6,'LineWidth',1.2,'HandleVisibility','off');
    text(convIter_p2,fitnessHistory_p2(convIter_p2),sprintf('  i%d\n  %.1fs',convIter_p2,timeToConv_p2),'FontSize',6,'Color',cPass2,'FontName','Times New Roman');
end
xlabel('Iteration','FontSize',8); ylabel('Mean Fitness (RMSE)','FontSize',8);
title('Convergence Curve','FontWeight','normal','FontSize',9);
legend('Location','northeast','FontSize',7);
set(gca,'FontName','Times New Roman','FontSize',8); grid on;

subplot(1,2,2);
cumT1=cumsum(timePerIter_p1); cumT2=cumsum(timePerIter_p2);
yyaxis left;
b_t1=bar(iP1-0.15,timePerIter_p1,0.3,'FaceColor',cPass1,'FaceAlpha',0.65,'EdgeColor',cPass1); hold on;
b_t2=bar(iP2+0.15,timePerIter_p2,0.3,'FaceColor',cPass2,'FaceAlpha',0.65,'EdgeColor',cPass2);
ylabel('Time per Iteration (s)','FontSize',8);
yyaxis right;
plot(iP1,cumT1,'-o','Color',cPass1*0.65,'LineWidth',0.9,'MarkerSize',4,'DisplayName','Single cumul');
plot(iP2,cumT2,'--s','Color',cPass2*0.65,'LineWidth',0.9,'MarkerSize',4,'DisplayName','Multi cumul');
ylabel('Cumulative Time (s)','FontSize',8);
xlabel('Iteration','FontSize',8);
title('Per-Iteration Timing','FontWeight','normal','FontSize',9);
legend([b_t1 b_t2],{'Single iter','Multi iter'},'Location','northwest','FontSize',6);
set(gca,'FontName','Times New Roman','FontSize',8); grid on;
sgtitle(sprintf('Convergence & Timing -- Case %s',caseDisplay),'FontWeight','normal','FontSize',10,'FontName','Times New Roman');
saveFig(fig3,[outBase '_Fitness.png'],[outBase '_Fitness.eps']);
fprintf('  Saved: _Fitness.png (convergence + per-iter timing)\n');

% Figure 4: Overall Performance Metrics
fig4=figure('Units','inches','Position',[0.5 0.5 7.16 4.5],'Color','w');
metricLabels={'Accuracy (MAE)','RMSE','Stability','Conv. Speed','R^2'};
barData=[agg.P1.MAE,agg.P1.RMSE,agg.P1.Stability,agg.P1.ConvSpeed,agg.P1.R2;...
         agg.P2.MAE,agg.P2.RMSE,agg.P2.Stability,agg.P2.ConvSpeed,agg.P2.R2];
nMods=size(barData,1); grpW=0.75; barW=grpW/nMods;
offs=linspace(-(grpW-barW)/2,(grpW-barW)/2,nMods);
for m=1:nMods
    xPos=(1:5)+offs(m);
    bar(xPos,barData(m,:),barW,'FaceColor',barColors(m,:),'EdgeColor','k','LineWidth',0.6); hold on;
    for k=1:5
        text(xPos(k),barData(m,k)+0.01*max(abs(barData(:,k))+eps),sprintf('%.3f',barData(m,k)),'HorizontalAlignment','center','FontSize',5.5,'FontName','Times New Roman');
    end
end
set(gca,'XTick',1:5,'XTickLabel',metricLabels,'FontName','Times New Roman','FontSize',8); xtickangle(18);
ylabel('Metric Value','FontSize',9);
title(sprintf('Performance Metrics -- Case %s',caseDisplay),'FontWeight','normal','FontSize',10,'FontName','Times New Roman');
legH=gobjects(nMods,1);
for m=1:nMods, legH(m)=bar(NaN,NaN,'FaceColor',barColors(m,:),'EdgeColor','k','LineWidth',0.6); end
legend(legH,modelLabels,'Location','northeast','FontSize',7); grid on; box on;
saveFig(fig4,[outBase '_Metrics.png'],[outBase '_Metrics.eps']);
fprintf('  Saved: _Metrics.png\n');

% Figure 5: Per-State RMSE and R2
fig5=figure('Units','inches','Position',[0.5 0.5 7.16 3.5],'Color','w');
subplot(1,2,1);
psRMSE=[metricsP1(:,1),metricsP2(:,1)];
b1=bar(psRMSE,'grouped');
for m=1:size(psRMSE,2), b1(m).FaceColor=barColors(m,:); b1(m).EdgeColor='k'; b1(m).LineWidth=0.5; end
set(gca,'XTickLabel',stateNames,'FontName','Times New Roman','FontSize',8); xtickangle(15);
ylabel('RMSE','FontSize',8); title('Per-State RMSE','FontWeight','normal','FontSize',9);
legend(modelLabels,'Location','best','FontSize',7); grid on;
subplot(1,2,2);
psR2=[metricsP1(:,3),metricsP2(:,3)];
b2=bar(psR2,'grouped');
for m=1:size(psR2,2), b2(m).FaceColor=barColors(m,:); b2(m).EdgeColor='k'; b2(m).LineWidth=0.5; end
set(gca,'XTickLabel',stateNames,'FontName','Times New Roman','FontSize',8); xtickangle(15);
ylabel('R^2','FontSize',8); title('Per-State R^2','FontWeight','normal','FontSize',9);
legend(modelLabels,'Location','best','FontSize',7); yline(0,'k--','LineWidth',0.5); grid on;
sgtitle(sprintf('Per-State Metrics -- Case %s',caseDisplay),'FontWeight','normal','FontSize',10,'FontName','Times New Roman');
saveFig(fig5,[outBase '_PerState.png'],[outBase '_PerState.eps']);
fprintf('  Saved: _PerState.png\n');

%% Excel Performance Spreadsheets
xlsFile=[outBase '_metrics.xlsx'];
if exist(xlsFile,'file'), delete(xlsFile); end

% Sheet 1: Aggregate Metrics with Timing
hdrAgg={'Model','MAE','RMSE','Stability','ConvFrac','R2','TotalTime_s','MeanIterTime_s','MinIterTime_s','MaxIterTime_s','TimeToConv_s','dMAE_vs_Single','dRMSE_vs_Single','dR2_vs_Single','Improvement'};
rowsAgg={};
mdefs={'Single-Step GWO (Pass-1)',agg.P1.MAE,agg.P1.RMSE,agg.P1.Stability,agg.P1.ConvSpeed,agg.P1.R2,agg.P1.TotalTime_s,agg.P1.MeanIterTime,agg.P1.MinIterTime,agg.P1.MaxIterTime,agg.P1.TimeToConv_s;
       'Multi-Step GWO (Pass-2)',agg.P2.MAE,agg.P2.RMSE,agg.P2.Stability,agg.P2.ConvSpeed,agg.P2.R2,agg.P2.TotalTime_s,agg.P2.MeanIterTime,agg.P2.MinIterTime,agg.P2.MaxIterTime,agg.P2.TimeToConv_s};

dM=mdefs{2,2}-mdefs{1,2}; dR=mdefs{2,3}-mdefs{1,3}; dR2=mdefs{2,6}-mdefs{1,6};
if dR2>0.05, verd='SIGNIFICANT IMPROVEMENT';
elseif dR2>0.01, verd='MODERATE IMPROVEMENT';
elseif dR2>-0.01, verd='COMPARABLE';
else, verd='NO IMPROVEMENT'; end

rowsAgg(1,:)={mdefs{1,1},mdefs{1,2},mdefs{1,3},mdefs{1,4},mdefs{1,5},mdefs{1,6},...
    mdefs{1,7},mdefs{1,8},mdefs{1,9},mdefs{1,10},mdefs{1,11},0,0,0,'Baseline'};
rowsAgg(2,:)={mdefs{2,1},mdefs{2,2},mdefs{2,3},mdefs{2,4},mdefs{2,5},mdefs{2,6},...
    mdefs{2,7},mdefs{2,8},mdefs{2,9},mdefs{2,10},mdefs{2,11},dM,dR,dR2,verd};

writetable(cell2table(rowsAgg,'VariableNames',hdrAgg),xlsFile,'Sheet','Aggregate_Metrics');

% Sheet 2: Per-State Breakdown
hdrPS={'Model','State','RMSE','MAE','R2','dRMSE_vs_Single','dR2_vs_Single'};
rowsPS={};
for s=1:nStates
    sn=stateNames{s};
    rowsPS(end+1,:)={'Single-Step GWO (Pass-1)',sn,metricsP1(s,1),metricsP1(s,2),metricsP1(s,3),0,0}; %#ok<AGROW>
    rowsPS(end+1,:)={'Multi-Step GWO (Pass-2)',sn,metricsP2(s,1),metricsP2(s,2),metricsP2(s,3),metricsP2(s,1)-metricsP1(s,1),metricsP2(s,3)-metricsP1(s,3)}; %#ok<AGROW>
end
writetable(cell2table(rowsPS,'VariableNames',hdrPS),xlsFile,'Sheet','PerState_Detail');

% Sheet 3: Convergence History
padN=MAX_ITER;
fH1=padFitnessVec(fitnessHistory_p1,padN); fH2=padFitnessVec(fitnessHistory_p2,padN);
iterCols=arrayfun(@(i)sprintf('Iter%d_Fit',i),1:padN,'UniformOutput',false);
hdrConv=[{'Model','NumIter','ConvIter','ConvFrac','FinalFit'},iterCols];
rowsConv={};
rowsConv(end+1,:)=[{'Single-Step GWO (Pass-1)',numel(fitnessHistory_p1),convIter_p1,convFrac_p1,fitnessHistory_p1(end)},num2cell(fH1)];
rowsConv(end+1,:)=[{'Multi-Step GWO (Pass-2)',numel(fitnessHistory_p2),convIter_p2,convFrac_p2,fitnessHistory_p2(end)},num2cell(fH2)];
writetable(cell2table(rowsConv,'VariableNames',hdrConv),xlsFile,'Sheet','Convergence_Fitness');

% Sheet 4: Execution Timing Breakdown
tH1=padFitnessVec(timePerIter_p1,padN); tH2=padFitnessVec(timePerIter_p2,padN);
timeCols=arrayfun(@(i)sprintf('Iter%d_Time_s',i),1:padN,'UniformOutput',false);
hdrT=[{'Model','TotalTime_s','MeanIter_s','MinIter_s','MaxIter_s','TimeToConv_s','ConvIter','ConvFrac'},timeCols];
rowsT={};
rowsT(end+1,:)=[{'Single-Step GWO (Pass-1)',agg.P1.TotalTime_s,agg.P1.MeanIterTime,agg.P1.MinIterTime,agg.P1.MaxIterTime,agg.P1.TimeToConv_s,convIter_p1,convFrac_p1},num2cell(tH1)];
rowsT(end+1,:)=[{'Multi-Step GWO (Pass-2)',agg.P2.TotalTime_s,agg.P2.MeanIterTime,agg.P2.MinIterTime,agg.P2.MaxIterTime,agg.P2.TimeToConv_s,convIter_p2,convFrac_p2},num2cell(tH2)];
writetable(cell2table(rowsT,'VariableNames',hdrT),xlsFile,'Sheet','Timing_Detail');

% Sheet 5: Bounds Comparison
hdrBnd={'Param','LB_orig','UB_orig','LB_refined','UB_refined','BestP1','BestP2','Delta'};
rowsBnd={};
for k=1:PARAM_COUNT
    rowsBnd(end+1,:)={sprintf('p%d',k),lb_orig(k),ub_orig(k),lb_p2(k),ub_p2(k),globalBest_p1(k),globalBest_p2(k),globalBest_p2(k)-globalBest_p1(k)}; %#ok<AGROW>
end
writetable(cell2table(rowsBnd,'VariableNames',hdrBnd),xlsFile,'Sheet','SearchBounds');
fprintf('\n  Excel saved: GWO_MultiStep_Case23_metrics.xlsx (5 sheets)\n');

%% Save Workspace
save([outBase '_workspace.mat'],'globalBest_p1','globalBest_p2',...
    'fitnessHistory_p1','fitnessHistory_p2',...
    'timePerIter_p1','timePerIter_p2','elapsedP1','elapsedP2',...
    'timeToConv_p1','timeToConv_p2',...
    'lb_orig','ub_orig','lb_p2','ub_p2',...
    'X_real','X_sim_p1','X_sim_p2','time',...
    'metricsP1','metricsP2','agg',...
    'convIter_p1','convIter_p2','convFrac_p1','convFrac_p2',...
    'CASE_NUM','SHRINK_FACTOR','POP_SIZE','MAX_ITER');
fprintf('  Workspace saved.\n\n');
fprintf('================================================================\n');
fprintf('  DONE -- outputs in: output/\n');
fprintf('================================================================\n\n');
diary off;
end

%% Core GWO Optimization Loop
function [globalBest,fitnessHist,convIter,timePerIter]=runGWO(staticData,lb,ub,popSize,maxIter,paramCount)
    wolves=InitializePopulation(popSize,paramCount,lb,ub);
    fitnessHist=zeros(maxIter,1);
    timePerIter=zeros(maxIter,1);
    convIter=maxIter;
    alphaWolf=wolves(1,:);
    for iter=1:maxIter
        t_it=tic;
        fitnessVals=EvaluateFitness24(wolves,staticData);
        [fitnessVals,idx]=sort(fitnessVals);
        wolves=wolves(idx,:);
        alphaWolf=wolves(1,:); betaWolf=wolves(2,:); deltaWolf=wolves(3,:);
        a=2-iter*(2/maxIter);
        for i=1:popSize
            for j=1:paramCount
                r1=rand;r2=rand;A1=2*a*r1-a;C1=2*r2;
                X1=alphaWolf(j)-A1*abs(C1*alphaWolf(j)-wolves(i,j));
                r1=rand;r2=rand;A2=2*a*r1-a;C2=2*r2;
                X2=betaWolf(j)-A2*abs(C2*betaWolf(j)-wolves(i,j));
                r1=rand;r2=rand;A3=2*a*r1-a;C3=2*r2;
                X3=deltaWolf(j)-A3*abs(C3*deltaWolf(j)-wolves(i,j));
                wolves(i,j)=max(lb(j),min(ub(j),(X1+X2+X3)/3));
            end
        end
        fitnessHist(iter)=mean(fitnessVals);
        timePerIter(iter)=toc(t_it);
        fprintf('    Iter %2d/%d | MeanFit=%.6f | Time=%.2fs\n',iter,maxIter,fitnessHist(iter),timePerIter(iter));
        if iter>1 && abs(fitnessHist(iter)-fitnessHist(iter-1))<1e-4
            fitnessHist=fitnessHist(1:iter);
            timePerIter=timePerIter(1:iter);
            fprintf('    Early stop iter=%d | total=%.2fs\n',iter,sum(timePerIter));
            break;
        end
    end
    globalBest=alphaWolf;
    timePerIter=timePerIter(1:numel(fitnessHist));
    finalF=fitnessHist(end); initF=fitnessHist(1);
    thr=finalF+0.01*abs(initF-finalF);
    convIter=numel(fitnessHist);
    for it=1:numel(fitnessHist)
        if fitnessHist(it)<=thr, convIter=it; break; end
    end
end

%% Longitudinal ODE Numerical Integration
function X_sim=simulateODE(params,uo_mean,time,u_data,de_data,fp_data,fn_data,fd_data,thrust,X0)
    [A_id,B_id]=formatParameters24(params,uo_mean);
    odefun=@(t,x)LongODE_tv(t,x,A_id,B_id,time,u_data,de_data,fp_data,fn_data,fd_data,thrust);
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t_sim,X_raw]=ode15s(odefun,time,X0,opts);
    X_sim=interp1(t_sim,X_raw,time,'linear','extrap');
end

function dx=LongODE_tv(t,x,A,B,time,u,elev,fp,fn,fd,thrust)
    de=interp1(time,elev,t,'linear','extrap');
    fpp=interp1(time,fp,t,'linear','extrap');
    fnn=interp1(time,fn,t,'linear','extrap');
    fdd=interp1(time,fd,t,'linear','extrap');
    ut=interp1(time,u,t,'linear','extrap');
    U=[de;thrust;fpp;fnn;fdd]; Atv=A; Atv(5,4)=ut; dx=Atv*x+B*U;
end

%% Performance Metrics
function m=computeMetrics3(y,yhat)
    y=y(:); yhat=yhat(:); ok=isfinite(y)&isfinite(yhat); y=y(ok); yhat=yhat(ok);
    if isempty(y), m=[NaN NaN NaN]; return; end
    err=y-yhat; rmse=sqrt(mean(err.^2)); mae=mean(abs(err));
    SS_res=sum(err.^2); SS_tot=sum((y-mean(y)).^2);
    r2=NaN; if SS_tot>0, r2=1-SS_res/SS_tot; end
    m=[rmse,mae,r2];
end

function v=padFitnessVec(fh,n)
    v=NaN(1,n); v(1:numel(fh))=fh(:)';
end

%% IEEE Figure Configuration
function applyIEEEDefaults()
    set(0,'DefaultFigureVisible','off','DefaultFigureColor','w',...
        'DefaultAxesColor','w','DefaultAxesFontName','Times New Roman',...
        'DefaultAxesFontSize',9,'DefaultAxesLineWidth',0.8,...
        'DefaultAxesTickDir','in','DefaultAxesBox','on',...
        'DefaultLineLineWidth',0.8,'DefaultTextFontName','Times New Roman',...
        'DefaultTextFontSize',9,'DefaultLegendFontName','Times New Roman',...
        'DefaultLegendFontSize',7,'DefaultTextColor','k',...
        'DefaultAxesXColor','k','DefaultAxesYColor','k');
end

function saveFig(fig,pngPath,epsPath)
    try, exportgraphics(fig,pngPath,'Resolution',600); catch, print(fig,pngPath,'-dpng','-r600'); end
    try, exportgraphics(fig,epsPath,'ContentType','vector'); catch, print(fig,epsPath,'-depsc','-r600'); end
    close(fig);
end
