function [AMatrix, BMatrix] = formatParameters6DOF(parameters, uo)
% formatParameters6DOF  Assembles continuous-time 12x12 A and 12x7 B matrices.
%
% State vector:
%   x = [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]^T
%
% Control input vector:
%   u = [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]^T
%
% Parameter vector (60 elements):
%   p(1:18)  - Stability derivatives (Xu,Xw,Xq, Yv,Yp,Yr, Zu,Zw,Zq, Lv,Lp,Lr, Mu,Mw,Mq, Nv,Np,Nr)
%   p(19:60) - Control derivatives (6 aerodynamic rows x 7 control inputs)

g = 9.81;

% Stability derivatives
Xu = parameters(1);  Xw = parameters(2);  Xq = parameters(3);
Yv = parameters(4);  Yp = parameters(5);  Yr = parameters(6);
Zu = parameters(7);  Zw = parameters(8);  Zq = parameters(9);
Lv = parameters(10); Lp = parameters(11); Lr = parameters(12);
Mu = parameters(13); Mw = parameters(14); Mq = parameters(15);
Nv = parameters(16); Np = parameters(17); Nr = parameters(18);

%% State Matrix A (12x12)
% States: [u, v, w, p, q, r, phi, theta, psi, xe, ye, h]
AMatrix = zeros(12, 12);

% Row 1: u_dot = Xu*u + Xw*w + Xq*q - g*theta (axial force dynamics)
AMatrix(1,  :) = [Xu,  0,  Xw,  0,  Xq,  0,   0,  -g,  0,  0,  0,  0];

% Row 2: v_dot = Yv*v + Yp*p + (Yr - uo)*r + g*phi (side force & kinematic coupling)
AMatrix(2,  :) = [0,  Yv,   0,  Yp,  0,  Yr - uo,  g,   0,  0,  0,  0,  0];

% Row 3: w_dot = Zu*u + Zw*w + (Zq + uo)*q (normal force & pitch coupling)
AMatrix(3,  :) = [Zu,  0,  Zw,  0,  Zq + uo,  0,   0,   0,  0,  0,  0,  0];

% Row 4: p_dot = Lv*v + Lp*p + Lr*r (rolling moment dynamics)
AMatrix(4,  :) = [0,  Lv,   0,  Lp,  0,  Lr,  0,   0,  0,  0,  0,  0];

% Row 5: q_dot = Mu*u + Mw*w + Mq*q (pitching moment dynamics)
AMatrix(5,  :) = [Mu,  0,  Mw,  0,  Mq,  0,   0,   0,  0,  0,  0,  0];

% Row 6: r_dot = Nv*v + Np*p + Nr*r (yawing moment dynamics)
AMatrix(6,  :) = [0,  Nv,   0,  Np,  0,  Nr,  0,   0,  0,  0,  0,  0];

% Row 7: phi_dot = p (linearized roll kinematics)
AMatrix(7,  :) = [0,   0,   0,   1,  0,   0,   0,   0,  0,  0,  0,  0];

% Row 8: theta_dot = q (pitch kinematics)
AMatrix(8,  :) = [0,   0,   0,   0,  1,   0,   0,   0,  0,  0,  0,  0];

% Row 9: psi_dot = r (linearized yaw kinematics)
AMatrix(9,  :) = [0,   0,   0,   0,  0,   1,   0,   0,  0,  0,  0,  0];

% Row 10: xe_dot = u (Earth-fixed North displacement)
AMatrix(10, :) = [1,   0,   0,   0,  0,   0,   0,   0,  0,  0,  0,  0];

% Row 11: ye_dot = v (Earth-fixed East displacement)
AMatrix(11, :) = [0,   1,   0,   0,  0,   0,   0,   0,  0,  0,  0,  0];

% Row 12: h_dot = -w + uo*theta (linearized climb rate)
AMatrix(12, :) = [0,   0,  -1,   0,  0,   0,   0,  uo,  0,  0,  0,  0];

%% Control Matrix B (12x7)
% Inputs: [delta_e, thrust, flapPos, flapNeg, flapDiff, aileron, rudder]
BMatrix = zeros(12, 7);
BMatrix(1,  :) = parameters(19:25);   % u_dot control row
BMatrix(2,  :) = parameters(26:32);   % v_dot control row
BMatrix(3,  :) = parameters(33:39);   % w_dot control row
BMatrix(4,  :) = parameters(40:46);   % p_dot control row
BMatrix(5,  :) = parameters(47:53);   % q_dot control row
BMatrix(6,  :) = parameters(54:60);   % r_dot control row
% Rows 7-12: pure kinematics (no direct actuator feedthrough)

end
