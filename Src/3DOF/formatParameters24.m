function [AMatrix,BMatrix] = formatParameters24(parameters,uo)
% Assembles 5-state longitudinal plant (A) and multi-control (B) matrices
% for states x = [u; w; q; theta; h] and inputs u_ctrl = [delta_e; delta_T; flapPos; flapNeg; flapDiff].

g = 9.81;

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

% 5-state longitudinal plant matrix
AMatrix = [Xn Xw Xq -g 0; Zu Zw Zq 0 0; Mu Mw Mq 0 0; 0 0 1 0 0; 0 -1 0 uo 0 ];

% Primary control derivatives (elevator and thrust)
X_delta_e = parameters(10);
X_delta_t = parameters(11);
Z_delta_e = parameters(12);
Z_delta_t = parameters(13);
M_delta_e = parameters(14);
M_delta_t = parameters(15);

% Auxiliary flap control derivatives
X_delta_flapPos = parameters(16);
X_delta_flapNeg = parameters(17);
X_delta_flapDiff = parameters(18);
Z_delta_flapPos = parameters(19);
Z_delta_flapNeg = parameters(20);
Z_delta_flapDiff = parameters(21);
M_delta_flapPos = parameters(22);
M_delta_flapNeg = parameters(23);
M_delta_flapDiff = parameters(24);

% 5x5 control allocation matrix
BMatrix = [X_delta_e X_delta_t X_delta_flapPos X_delta_flapNeg X_delta_flapDiff ;...
    Z_delta_e Z_delta_t Z_delta_flapPos Z_delta_flapNeg Z_delta_flapDiff ;...
    M_delta_e M_delta_t M_delta_flapPos M_delta_flapNeg M_delta_flapDiff ;...
    0 0 0 0 0 ;
    0 0 0 0 0 ];

end
