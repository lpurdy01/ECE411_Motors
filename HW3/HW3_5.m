%% Problem 5 stuff
% Finds: torque angle δ; Eq, Ed; stator current I0; PF angle θ; Iq, Id.
% Assumptions: V on q-axis (v_d = 0), r_s = 0, uniform L_s.

clear; clc;

%% Givens (from problem)
f_e   = 60;                 % electrical frequency [Hz]
we    = 2*pi*f_e;           % electrical speed [rad/s]
Ls    = 26e-3;              % synchronous inductance [H]
P_mech= 20e3;               % shaft power (motoring) [W]
V_RMS = 277;                % stator line-neutral RMS voltage [V_LN_RMS]
E_RMS = 310;                % back-EMF line-neutral RMS [V_LN_RMS]

%% Convert RMS -> peak for dq (class notes: qd values correspond to LN peak)
V0 = sqrt(2)*V_RMS;         % peak line-neutral stator voltage
E0 = sqrt(2)*E_RMS;         % peak line-neutral back-EMF

%% Solve for torque angle δ using power balance:
% P = (3/2) * V0 * Iq  and  Iq = E0*sin(δ)/(we*Ls)
% => sin(δ) = (2*P*we*Ls)/(3*V0*E0)
sin_delta = (2*P_mech*we*Ls)/(3*V0*E0);
sin_delta = max(-1,min(1,sin_delta));   % clamp for numerical safety
delta = asin(sin_delta);                % rad, choose principal (motoring 0<δ<π/2)

%% Back-EMF components (Lipo–Novotny)
Eq =  E0*cos(delta);
Ed = -E0*sin(delta);

%% Current components from V = E + j*we*Ls*I   (vd=0, vq=V0)
Iq = -Ed/(we*Ls);                       % =  E0*sin(δ)/(we*Ls)
Id = (Eq - V0)/(we*Ls);                 % = (E0*cos(δ) - V0)/(we*Ls)

%% Magnitude & PF angle (LN convention: θ = atan2(-Id, Iq))
I0 = hypot(Iq, Id);                  % peak line-neutral current
theta = atan2(-Id, Iq);                 % rad, positive for lagging

%% Print results
fprintf('--- Problem 5a results (Lipo–Novotny, dq peaks) ---\n');
fprintf('δ (torque angle)  = %.6f rad  (%.3f deg)\n', delta, rad2deg(delta));
fprintf('Eq          = %.3f V\n', Eq);
fprintf('Ed         = %.3f V\n', Ed);
fprintf('Iq          = %.3f A\n', Iq);
fprintf('Id          = %.3f A\n', Id);
fprintf('I0          = %.3f A\n', I0);
fprintf('θ (PF angle)      = %.6f rad  (%.3f deg)\n', theta, rad2deg(theta));

%% Optional quick checks
% P_check = (3/2)*V0*Iq;  % should be ~P_mech
% fprintf('Power check (from dq): %.2f W\n', P_check);
