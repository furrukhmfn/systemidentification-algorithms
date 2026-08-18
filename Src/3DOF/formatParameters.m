function [AMatrix,BMatrix] = formatParameters(parameters,uo, valueOfG)
% Constructs linearized 5-state longitudinal state (A) and control (B) matrices
% for states x = [u; w; q; theta; h] and controls u_ctrl = [delta_e; delta_T].

g = valueOfG;

% Dimensional stability derivatives
Xn = parameters(1);  % Xu
Xw = parameters(2);
Xq = parameters(3);
Zu = parameters(4);
Zw = parameters(5);
Zq = parameters(6);
Mu = parameters(7);
Mw = parameters(8);
Mq = parameters(9);

% 5-state longitudinal system matrix
AMatrix = [Xn Xw Xq -g 0; Zu Zw Zq 0 0; Mu Mw Mq 0 0; 0 0 1 0 0; 0 -1 0 uo 0 ];

% Control effectiveness derivatives
X_delta_e = parameters(10);
X_delta_t = parameters(11);
Z_delta_e = parameters(12);
Z_delta_t = parameters(13);
M_delta_e = parameters(14);
M_delta_t = parameters(15);

% Control input matrix [delta_e, delta_T]
BMatrix = [X_delta_e X_delta_t; Z_delta_e Z_delta_t; M_delta_e M_delta_t; 0 0; 0 0];

end