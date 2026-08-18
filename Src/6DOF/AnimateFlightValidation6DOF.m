function AnimateFlightValidation6DOF(outputFileName, caseNumber, algorithmName, varargin)
% AnimateFlightValidation6DOF  3D trajectory animation comparing flight data vs identified model.
% Integrates identified 6-DOF state-space equations via ODE15s and renders simultaneous 3D
% aircraft silhouettes (measured vs estimated) alongside live telemetry panels.
%
% Usage:
%   AnimateFlightValidation6DOF('output/GWO_Case41_6DOF.mat', 41, 'GWO')
%   AnimateFlightValidation6DOF(..., 'SaveVideo', true, 'Subsample', 8)

%% Options and Parsing
p = inputParser;
addParameter(p, 'SaveVideo', true);
addParameter(p, 'FrameRate', 20);
addParameter(p, 'Subsample', 8);
addParameter(p, 'TrailFrac', 0.22);
parse(p, varargin{:});
doSave    = p.Results.SaveVideo;
fps       = p.Results.FrameRate;
sub       = p.Results.Subsample;
trailFrac = p.Results.TrailFrac;

%% Load Identification Results
fprintf('Loading: %s\n', outputFileName);
ld           = load(outputFileName);
Vb           = ld.Vb;
pqr_raw      = ld.pqr;
pts          = ld.phi_theta_psi;
Xe           = ld.Xe;
globalBest   = ld.globalBest;
sd           = ld.staticData;
caseDisp     = sprintf('(%d,%d)', floor(caseNumber/10), mod(caseNumber,10));

%% Reference Flight Telemetry
try; mgn = ld.inputData.ValidationMargin; catch; mgn = 0; end
N  = size(Vb,1);  dS = 1+mgn;  dE = N-mgn;

u_r   = Vb(dS:dE,1);
phi_r = pts(1,dS:dE)';
the_r = pts(2,dS:dE)';
psi_r = pts(3,dS:dE)';
xe_r  = Xe(dS:dE,1);
ye_r  = Xe(dS:dE,2);
h_r   = -Xe(dS:dE,3);

t_all = sd.time_all;
Nobs  = length(u_r);

% 12-state flight telemetry matrix
X_real = [u_r, sd.v_all, sd.w_all, sd.p_all, sd.q_all, sd.r_all, ...
          phi_r, the_r, psi_r, xe_r, ye_r, h_r]';

%% State-Space Matrices
u0   = sd.u0_trim;
xTrm = sd.x_trim;
uTrm = sd.u_trim;
[A, B] = formatParameters6DOF(globalBest, u0);

%% Actuator Deflection Perturbations
try
    ede = deg2rad(ld.elevator(dS:dE,1) - uTrm(1));
    efp = deg2rad(ld.flapPos( dS:dE,1) - uTrm(3));
    efn = deg2rad(ld.flapNeg( dS:dE,1) - uTrm(4));
    efd = deg2rad(ld.flapDiff(dS:dE,1) - uTrm(5));
    eda = deg2rad(ld.aileron( dS:dE,1) - uTrm(6));
    edr = deg2rad(ld.rudder(  dS:dE,1) - uTrm(7));
    nIn = 7;
catch
    ede = deg2rad(ld.elevator(dS:dE,1) - uTrm(1));
    efp = deg2rad(ld.flapPos( dS:dE,1) - uTrm(3));
    efn = deg2rad(ld.flapNeg( dS:dE,1) - uTrm(4));
    efd = deg2rad(ld.flapDiff(dS:dE,1) - uTrm(5));
    eda = zeros(Nobs,1);
    edr = zeros(Nobs,1);
    nIn = 5;
end

%% Perturbation Model Integration (ODE15s)
fprintf('Integrating with ODE15s ...\n');
X0p = X_real(:,1) - xTrm;
odef = @(t,dx) afODE(t, dx, A, B, t_all, u_r, ede, efp, efn, efd, eda, edr, nIn);
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
try
    [ts, Xr] = ode15s(odef, t_all, X0p, opts);
    Xp = interp1(ts, Xr, t_all, 'linear', 'extrap');
    Xe2 = (Xp + xTrm')';
    fprintf('  ODE15s OK.\n');
catch ME
    warning('ODE15s failed (%s). Falling back to Euler integration.', ME.message);
    Xp = zeros(12, Nobs);  Xp(:,1) = X0p;
    dt = mean(diff(t_all));
    for k = 1:Nobs-1
        Atv = A;  Atv(12,8) = u_r(k);
        if nIn==7; U=[ede(k);0;efp(k);efn(k);efd(k);eda(k);edr(k)];
        else;      U=[ede(k);0;efp(k);efn(k);efd(k)]; end
        Xp(:,k+1) = Xp(:,k) + dt*(Atv*Xp(:,k) + B*U);
    end
    Xe2 = Xp + xTrm;
end

% Bounding envelope protection for display rendering
for row = 1:12
    rm = mean(X_real(row,:));
    rs = max(max(X_real(row,:)) - min(X_real(row,:)), 0.01);
    Xe2(row,:) = min(max(Xe2(row,:), rm - 4*rs), rm + 4*rs);
end

xe_e  = Xe2(10,:)';   ye_e  = Xe2(11,:)';   h_e   = Xe2(12,:)';
phi_e = Xe2(7, :)';   the_e = Xe2(8, :)';   psi_e = Xe2(9, :)';
u_e   = Xe2(1, :)';

%% Subsampling for Smooth Playback
idx = 1:sub:Nobs;
Nfr = length(idx);
sXR = xe_r(idx);  sYR = ye_r(idx);  sHR = h_r(idx);
sXE = xe_e(idx);  sYE = ye_e(idx);  sHE = h_e(idx);
sPR = phi_r(idx); sTR = the_r(idx); sWR = psi_r(idx);
sPE = phi_e(idx); sTE = the_e(idx); sWE = psi_e(idx);
sUR = u_r(idx);   sUE = u_e(idx);
sFT = t_all(idx);

%% Aircraft 3D Model Scaling
extY = max(ye_r) - min(ye_r);
extH = max(h_r)  - min(h_r);
extX = max(xe_r) - min(xe_r);
cRef = max(min(extY, extH), 20);
acS  = min(cRef * 0.12, 22);
acS  = max(acS, 2);
fprintf('Aircraft icon: %.1f m  (extY=%.0f m  extH=%.0f m  extX=%.0f m)\n', ...
    acS, extY, extH, extX);

%% Aircraft Silhouette Vertices (Body Frame: +X Nose, +Y Right, +Z Up)
% Fuselage and swept wing geometry
bX = acS * [ 1.30, 1.05, 0.70, 0.35, 0.15, ...
            -0.25, -0.58, -0.48, ...
            -0.80, -1.10, -1.20, ...
            -1.10, -0.80, ...
            -0.48, -0.58, -0.25, ...
             0.15,  0.35,  0.70, 1.05, 1.30];
bY = acS * [ 0.00,  0.05,  0.08,  0.10,  0.10, ...
             0.75,   0.75,  0.10, ...
             0.06,   0.04,  0.00, ...
            -0.04,  -0.06, ...
            -0.10,  -0.75, -0.75, ...
            -0.10,  -0.10, -0.08, -0.05,  0.00];
bZ = zeros(1,21);
bV = [bX(:), bY(:), bZ(:)];

% Horizontal and vertical stabilizers
rHX = acS * [-0.78, -1.05, -1.12, -1.10];
rHY = acS * [ 0.05,  0.05,  0.28,  0.27];
rHZ = zeros(1,4);
rHV = [rHX(:), rHY(:), rHZ(:)];

lHX = acS * [-0.78, -1.05, -1.12, -1.10];
lHY = acS * [-0.05, -0.05, -0.28, -0.27];
lHZ = zeros(1,4);
lHV = [lHX(:), lHY(:), lHZ(:)];

vSX = acS * [-0.78, -1.10, -1.00, -0.82];
vSY = zeros(1,4);
vSZ = acS * [ 0.00,  0.00,  0.28,  0.18];
vSV = [vSX(:), vSY(:), vSZ(:)];

%% Figure Layout & Graphics
BG  = [0.06, 0.06, 0.12];
PAN = [0.09, 0.09, 0.17];
AX  = [0.60, 0.65, 0.78];
GR  = [0.20, 0.20, 0.30];
CR  = [0.28, 0.62, 1.00];   % Flight truth (blue)
CE  = [1.00, 0.38, 0.28];   % Identified model (coral)
CW  = [0.97, 0.97, 0.97];

fig = figure('Color',BG,'Units','pixels','Position',[50 50 1280 720], ...
    'Name',sprintf('Flight Animation - %s %s',algorithmName,caseDisp), ...
    'NumberTitle','off','Resize','off');
set(fig,'InvertHardcopy','off');

% 3D Earth-Frame Viewport
ax3 = axes('Parent',fig,'Position',[0.02 0.08 0.57 0.88], ...
    'Color',PAN,'XColor',AX,'YColor',AX,'ZColor',AX, ...
    'GridColor',GR,'GridAlpha',0.45,'FontSize',8,'LineWidth',0.7);
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');

% Background full trajectory traces
plot3(ax3, xe_r, ye_r, h_r, '-','Color',[CR,0.38],'LineWidth',1.4);
plot3(ax3, xe_e, ye_e, h_e, '-','Color',[CE,0.38],'LineWidth',1.4);

% Animated flight trail lines
hTR = plot3(ax3, nan,nan,nan, '-','Color',CR,'LineWidth',3.0);
hTE = plot3(ax3, nan,nan,nan, '-','Color',CE,'LineWidth',3.0);

% Starting waypoints
plot3(ax3,xe_r(1),ye_r(1),h_r(1),'o','MarkerSize',8,...
    'MarkerFaceColor',CR,'MarkerEdgeColor',CW,'LineWidth',1.2);
plot3(ax3,xe_e(1),ye_e(1),h_e(1),'o','MarkerSize',8,...
    'MarkerFaceColor',CE,'MarkerEdgeColor',CW,'LineWidth',1.2);

xlabel(ax3,'X_e (m)','Color',AX,'FontSize',9);
ylabel(ax3,'Y_e (m)','Color',AX,'FontSize',9);
zlabel(ax3,'h (m)',  'Color',AX,'FontSize',9);
title(ax3, sprintf('%s  Case %s       Blue = Real       Red = Estimated', ...
    algorithmName, caseDisp), 'Color',CW,'FontSize',10,'FontWeight','normal');

view(ax3,-52,20); set(ax3,'Projection','perspective');
axis(ax3,'vis3d'); axis(ax3,'tight');

hTxt = text(ax3, xe_r(1), ye_r(1), h_r(1),'t = 0.0 s', ...
    'Color',CW,'FontSize',10,'FontWeight','bold','HorizontalAlignment','left');

[hBR,hLR,hRR,hVR] = afMkAC(ax3, bV, lHV, rHV, vSV, CR);
[hBE,hLE,hRE,hVE] = afMkAC(ax3, bV, lHV, rHV, vSV, CE);

% Live Telemetry Subpanels
pLbl = {'u  (m/s)', '\phi  (deg)', '\theta  (deg)', 'h  (m)'};
pR   = {u_r,          rad2deg(phi_r),   rad2deg(the_r),   h_r};
pE   = {u_e,          rad2deg(phi_e),   rad2deg(the_e),   h_e};
pfR  = {sUR,          rad2deg(sPR),     rad2deg(sTR),     sHR};
pfE  = {sUE,          rad2deg(sPE),     rad2deg(sTE),     sHE};
fT   = t_all(1:Nobs);

nP   = 4;
axP  = gobjects(nP,1);
hDR  = gobjects(nP,1);
hDE  = gobjects(nP,1);
hC   = gobjects(nP,1);
pH   = 0.180; pGap = 0.026; pX = 0.640; pW = 0.348;

for i = 1:nP
    yb = 0.055 + (nP-i)*(pH+pGap);
    axP(i) = axes('Parent',fig,'Position',[pX yb pW pH], ...
        'Color',PAN,'XColor',AX,'YColor',AX,'GridColor',GR, ...
        'GridAlpha',0.5,'FontSize',7.5,'LineWidth',0.6);
    hold(axP(i),'on'); grid(axP(i),'on'); box(axP(i),'on');

    plot(axP(i), fT, pR{i}, '-','Color',[CR,0.92],'LineWidth',1.5);
    axis(axP(i),'tight'); xlim(axP(i),[fT(1),fT(end)]);
    yl = ylim(axP(i)); pd = max((yl(2)-yl(1))*0.15, 1e-3);
    ylim(axP(i),[yl(1)-pd, yl(2)+pd]); yl = ylim(axP(i));
    set(axP(i),'YLimMode','manual');

    plot(axP(i), fT, pE{i}, '-','Color',[CE,0.88],'LineWidth',1.5);

    hC(i)  = plot(axP(i),[sFT(1),sFT(1)],yl,'w-','LineWidth',1.8);
    hDR(i) = plot(axP(i),sFT(1),pfR{i}(1),'o','MarkerSize',6, ...
        'MarkerFaceColor',CR,'MarkerEdgeColor',CW,'LineWidth',0.7);
    hDE(i) = plot(axP(i),sFT(1),pfE{i}(1),'o','MarkerSize',6, ...
        'MarkerFaceColor',CE,'MarkerEdgeColor',CW,'LineWidth',0.7);

    ylabel(axP(i),pLbl{i},'Color',AX,'FontSize',8,'Interpreter','tex');
    if i==nP; xlabel(axP(i),'Time (s)','Color',AX,'FontSize',8);
    else;     set(axP(i),'XTickLabel',{}); end
    xlim(axP(i),[fT(1),fT(end)]);
end
legend(axP(1),{'Real','Estimated'},'TextColor',CW,'Color',[0.10 0.10 0.18], ...
    'EdgeColor',GR,'FontSize',7,'Location','best');

%% Video Recording Setup
if doSave
    if ~exist('output','dir'); mkdir('output'); end
    vName = sprintf('output/FlightAnimation_%s_Case%d_6DOF.mp4',algorithmName,caseNumber);
    drawnow;
    vw = VideoWriter(vName,'MPEG-4');
    vw.FrameRate = fps;  vw.Quality = 95;
    open(vw);
    fprintf('Saving video: %s\n', vName);
end

%% Animation Playback Loop
tLen = max(3, round(Nfr * trailFrac));
fprintf('Rendering %d frames...\n', Nfr);

for fr = 1:Nfr
    f0 = max(1, fr-tLen);

    set(hTR,'XData',sXR(f0:fr),'YData',sYR(f0:fr),'ZData',sHR(f0:fr));
    set(hTE,'XData',sXE(f0:fr),'YData',sYE(f0:fr),'ZData',sHE(f0:fr));

    set(hTxt,'Position',[sXR(fr), sYR(fr), sHR(fr)+acS*2.2], ...
        'String', sprintf('t = %.1f s', sFT(fr)));

    % Rotate and translate measured aircraft icon
    Rr = afRot(sWR(fr), sTR(fr), sPR(fr));
    pR2 = [sXR(fr); sYR(fr); sHR(fr)];
    afMvAC(hBR,hLR,hRR,hVR, bV,lHV,rHV,vSV, Rr, pR2);

    % Rotate and translate estimated aircraft icon
    Re = afRot(sWE(fr), sTE(fr), sPE(fr));
    pE2 = [sXE(fr); sYE(fr); sHE(fr)];
    afMvAC(hBE,hLE,hRE,hVE, bV,lHV,rHV,vSV, Re, pE2);

    % Update panel tracking cursors
    for i = 1:nP
        yl = ylim(axP(i));
        set(hC(i), 'XData',[sFT(fr),sFT(fr)],'YData',yl);
        set(hDR(i),'XData',sFT(fr),'YData',pfR{i}(fr));
        set(hDE(i),'XData',sFT(fr),'YData',pfE{i}(fr));
    end

    drawnow limitrate;
    if doSave; writeVideo(vw, getframe(fig)); end
end

if doSave; close(vw); fprintf('Done: %s\n', vName); end

end

%% Local Helper Functions

function ddx = afODE(t, dx, A, B, tVec, u_abs, ede, efp, efn, efd, eda, edr, nIn)
% Continuous perturbation state-space derivative evaluation
de = interp1(tVec, ede, t,'linear','extrap');
fp = interp1(tVec, efp, t,'linear','extrap');
fn = interp1(tVec, efn, t,'linear','extrap');
fd = interp1(tVec, efd, t,'linear','extrap');
ut = interp1(tVec, u_abs, t,'linear','extrap');
if nIn == 7
    da = interp1(tVec, eda, t,'linear','extrap');
    dr = interp1(tVec, edr, t,'linear','extrap');
    U  = [de; 0; fp; fn; fd; da; dr];
else
    U  = [de; 0; fp; fn; fd];
end
Atv = A;  Atv(12,8) = ut;
ddx = Atv*dx + B*U;
end

function R = afRot(psi, theta, phi)
% Direction cosine matrix (Body to Earth coordinate transformation via ZYX Euler angles)
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta); cy=cos(psi); sy=sin(psi);
R = [cy*ct, cy*st*sp-sy*cp, cy*st*cp+sy*sp;
     sy*ct, sy*st*sp+cy*cp, sy*st*cp-cy*sp;
     -st,   ct*sp,          ct*cp];
end

function [hBody, hHL, hHR, hVS] = afMkAC(ax, bV, lHV, rHV, vSV, clr)
% Aircraft patch creation
hBody = patch(ax,'XData',bV(:,1), 'YData',bV(:,2), 'ZData',bV(:,3), ...
    'FaceColor',clr,      'EdgeColor',clr*0.55,'FaceAlpha',0.93,'LineWidth',0.6);
hHL   = patch(ax,'XData',lHV(:,1),'YData',lHV(:,2),'ZData',lHV(:,3), ...
    'FaceColor',clr*0.72, 'EdgeColor',clr*0.50,'FaceAlpha',0.90,'LineWidth',0.5);
hHR   = patch(ax,'XData',rHV(:,1),'YData',rHV(:,2),'ZData',rHV(:,3), ...
    'FaceColor',clr*0.72, 'EdgeColor',clr*0.50,'FaceAlpha',0.90,'LineWidth',0.5);
hVS   = patch(ax,'XData',vSV(:,1),'YData',vSV(:,2),'ZData',vSV(:,3), ...
    'FaceColor',clr*0.50, 'EdgeColor',clr*0.40,'FaceAlpha',0.88,'LineWidth',0.5);
end

function afMvAC(hB,hL,hR,hV, bV,lHV,rHV,vSV, Rot, pos)
% Transform all aircraft components
afSetP(hB, bV,  Rot, pos);
afSetP(hL, lHV, Rot, pos);
afSetP(hR, rHV, Rot, pos);
afSetP(hV, vSV, Rot, pos);
end

function afSetP(h, verts, R, pos)
% Apply 3D rotation and translation to vertex array
V = (R * verts')';
set(h,'XData',V(:,1)+pos(1),'YData',V(:,2)+pos(2),'ZData',V(:,3)+pos(3));
end
