% DC machine Parameters
% Based off a 25 hp, 500 rpm DC machine
K = 4.0;          % [V/(rad/s)] = [N*m/A] motor constant
R = 0.115;        % [Ohms] Hot resistance
L = 0.011;        % [H] DC machine inductance
J = 0.3;          % [kg*m^2] Moment of inertia
B = 1.0;          % [N*m/(rad/s)] Damping coefficient
w_step = 22.5;    % [rad/s] Open loop commanded speed
K_est = K;        % Estimated motor constant = actual
Tl = 0;
I_step = 5.5;
K_pa = 7;
K_ia = 220;

