%% problem 6 a
% given: Eq = E_o, Ed = 0; I leads E by gamma rs=0 pu; Xs=0.5 pu @ we=1; P=1 pu
% find: Vq, Vd, delta, I0, theta, Iq, Id

clear; clc;

%% givens (pu)
E_o = 1.20;   % internal back emf [pu]
gamma = 0.25;   % internal pf angle (I leads E) rad
Xs    = 0.5;    % sync react @ we=1 pu
P     = 1.0;    % rated power pu

%% current mag from power (LN scaling)
% P = (3/2)*E_inf*I0*cos(gamma)
I0 = (2/3)*(P/(E_o*cos(gamma)));

%% current comps (lipo-novotny: q=cos, d=-sin)
Iq =  I0*cos(gamma);
Id = -I0*sin(gamma);

%% stator voltage (rs=0): V = E + j*Xs*I
% Vq = -Xs*Id + Eq ; Vd = Xs*Iq + Ed, with Eq=E_inf, Ed=0
Vq = -Xs*Id + E_o;
Vd =  Xs*Iq;

%% torque angle delta = angle(V) from q-axis
delta = atan2(-Vd, Vq);  % q = real, d = -imag

%% pf angle theta = angle(I) - angle(V)
% angle(I) = atan2(-Id, Iq) = gamma
theta = gamma - delta;

%% prints
fprintf('p 6a (perunit)\n');
fprintf('delta = %.6f rad (%.3f deg)\n', delta, rad2deg(delta));
fprintf('Vq    = %.6f pu\n', Vq);
fprintf('Vd    = %.6f pu\n', Vd);
fprintf('I0    = %.6f pu\n', I0);
fprintf('theta = %.6f rad (%.3f deg)\n', theta, rad2deg(theta));
fprintf('Iq    = %.6f pu\n', Iq);
fprintf('Id    = %.6f pu\n', Id);

%% part b
% given: same E_inf and I0 as part (a)
% find gamma for max motoring P, then Pmax, Vq, Vd, delta, theta, Iq, Id


%% optimal gamma (max cos)
gamma_opt = 0;     % rad, current aligned with E

%% currents (lipo novotny map)
Iq_b =  I0*cos(gamma_opt);
Id_b = -I0*sin(gamma_opt);   % =0 at opt

%% terminal volt comps
Vq_b = -Xs*Id_b + E_o;     % = E_inf
Vd_b =  Xs*Iq_b;             % = Xs*I0

%% torque angle (angle of V)
delta_b = atan2(-Vd_b, Vq_b);

%% pf angle
theta_b = gamma_opt - delta_b;

%% power @ opt
Pmax = (3/2)*E_o*I0*cos(gamma_opt);

%% print stuff
fprintf('\n--- p 6b (max power fixed E_inf & I0) ---\n');
fprintf('gamma* = %.6f rad (%.3f deg)\n', gamma_opt, rad2deg(gamma_opt));
fprintf('Pmax   = %.6f pu\n', Pmax);
fprintf('Vq     = %.6f pu\n', Vq_b);
fprintf('Vd     = %.6f pu\n', Vd_b);
fprintf('delta  = %.6f rad (%.3f deg)\n', delta_b, rad2deg(delta_b));
fprintf('theta  = %.6f rad (%.3f deg)\n', theta_b, rad2deg(theta_b));
fprintf('Iq     = %.6f pu\n', Iq_b);
fprintf('Id     = %.6f pu\n', Id_b);