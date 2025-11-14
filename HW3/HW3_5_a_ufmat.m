%% problem 5 a
% find: delta, Eq/Ed, I0, PF angle, Iq/Id
% assume: V on q-axis (vd=0), rs=0, Ls uniform, lipo-novotny

clear; clc;

%% givens
f_e   = 60;              % electrical freq Hz
we    = 2*pi*f_e;        % elec speed rad/s
Ls    = 26e-3;           % sync inductance H
P_mech= 20e3;            % mech power W (motoring)
V_RMS = 277;             % line-neutral rms V
E_RMS = 310;             % back-emf ln rms V

%% rms -> peak (dq uses LN peak)
V0 = sqrt(2)*V_RMS;
E0 = sqrt(2)*E_RMS;

%% torque angle delta from power balance
% sin(delta) = (2*P*we*Ls)/(3*V0*E0)
sin_delta = (2*P_mech*we*Ls)/(3*V0*E0);
sin_delta = max(-1,min(1,sin_delta));
delta = asin(sin_delta);    % principal (0..pi/2 for motoring)

%% back-emf comps (lipo-novotny)
Eq =  E0*cos(delta);
Ed = -E0*sin(delta);

%% currents from V = E + j*we*Ls*I  vq=V0 vd=0
Iq = -Ed/(we*Ls);
Id = (Eq - V0)/(we*Ls);

%% magnitude + pf angle (ln conv: theta = atan2(-Id, Iq))
I0    = hypot(Iq, Id);      % peak line-neutral current
theta = atan2(-Id, Iq);     % + for lagging

%% print
fprintf('p 5a \n');
fprintf('delta = %.6f rad (%.3f deg)\n', delta, rad2deg(delta));
fprintf('Eq    = %.3f V\n', Eq);
fprintf('Ed    = %.3f V\n', Ed);
fprintf('Iq    = %.3f A\n', Iq);
fprintf('Id    = %.3f A\n', Id);
fprintf('I0    = %.3f A\n', I0);
fprintf('theta = %.6f rad (%.3f deg)\n', theta, rad2deg(theta));


