%% problem 5 b
% find: delta, Eq/Ed, I0, PF angle, Iq/Id
% assume: PF = 1  therefore  theta = 0,  Id = 0
%         V on q-axis (vd=0), rs=0, Ls uniform, lipo-novotny

clear; clc;

%% givens (same machine as part a)
f_e   = 60;              % electrical freq Hz
we    = 2*pi*f_e;        % elec speed rad/s
Ls    = 26e-3;           % sync inductance H
P_mech= 20e3;            % mech power W (motoring)
V_RMS = 277;             % terminal line-neutral rms V

%% rms -> peak (dq uses LN peak)
V0 = sqrt(2)*V_RMS;      % stator voltage peak (on q-axis)

%% PF = 1 => Id = 0 and current aligned with V (theta = 0)
theta = 0;               % rad
Id    = 0;

%% From power balance with PF=1:
% P = (3/2) * V0 * Iq   (since vd=0 and Id=0)
Iq = (2*P_mech) / (3*V0);

%% Back-emf components from stator eqn  V = E + j*we*Ls*I  
% v_d = we*Ls*Iq + E_d = 0  ->  E_d = -we*Ls*Iq
% v_q = -we*Ls*Id + E_q = V0 ->  E_q = V0
Ed = -we*Ls*Iq;
Eq =  V0;

%% Internal EMF magnitude and torque angle
E0    = hypot(Eq, Ed);
delta = atan2(-Ed, Eq);  % δ = angle from q-axis to E  (LN: Ed = -E0 sinδ)

%% Current magnitude peak line-neutral; RMS shown for reference)
I0 = hypot(Iq, Id);      % = |Iq| here

%% print
fprintf('p 5b (unity PF) \n');
fprintf('delta = %.6f rad (%.3f deg)\n', delta, rad2deg(delta));
fprintf('Eq    = %.3f V\n', Eq);
fprintf('Ed    = %.3f V\n', Ed);
fprintf('I0    = %.3f A  \n', I0);
fprintf('theta = %.6f rad (PF = 1)\n', theta);
fprintf('Iq    = %.3f A\n', Iq);
fprintf('Id    = %.3f A\n', Id);

%% check
P_chk = (3/2)*V0*Iq;
  fprintf('Power check = %.2f W\n', P_chk);
