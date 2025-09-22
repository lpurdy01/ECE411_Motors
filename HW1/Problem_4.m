% File: hw1_q4a_universal_PMDC_plots.m
% Purpose: Plot ideal steady-state PMDC motor behavior vs speed (ω) for two motors
% Plots: (1) Torque, (2) Armature current, (3) Armature copper loss, (4) Mechanical power
% Notes: Uses ideal 24 V battery, SI units, and clamps Ia>=0 (no regen in this ideal model)

clear; clc;

%% ---------------- Parameters ----------------
Vdc = 24;                    % [V] ideal battery

% Motor 1: C23_L55_20
mot(1).name = 'Motor 1: C23\_L55\_20';
mot(1).Ra   = 1.43;          % [ohm]
mot(1).K    = 0.0991;        % [Nm/A] == [V·s/rad]

% Motor 2: C34_L60_10
mot(2).name = 'Motor 2: C34\_L60\_10';
mot(2).Ra   = 0.43;          % [ohm]
mot(2).K    = 0.08;          % [Nm/A]

%% ------------- Speed axis (rad/s) -------------
% Use 0 → slightly above the larger no-load speed among the two motors.
w0 = Vdc ./ [mot.K];               % no-load speeds [rad/s] for each motor
w_max = 1.05 * max(w0);            % 5% margin beyond the highest no-load speed
w = linspace(0, w_max, 600).';     % column vector [rad/s]

%% ------------- Compute curves -------------
for i = 1:2
    K  = mot(i).K;
    Ra = mot(i).Ra;

    Ia = (Vdc - K.*w) ./ Ra;       % [A] ideal Ia vs speed
    Ia = max(Ia, 0);               % clamp negative current to 0 (no regen in this idealization)

    T  = K .* Ia;                   % [Nm] torque
    Ploss = Ia.^2 .* Ra;            % [W] armature copper loss
    Pmech = T .* w;                 % [W] mechanical output power

    % Store
    mot(i).Ia    = Ia;
    mot(i).T     = T;
    mot(i).Ploss = Ploss;
    mot(i).Pmech = Pmech;
end

%% ------------- Plotting -------------
figure('Name','4.a Ideal Motor and Battery Performance','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) Torque vs ω
nexttile;
plot(w, mot(1).T, 'LineWidth',1.8); hold on;
plot(w, mot(2).T, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Torque  T  [N·m]');
title('Torque vs. Speed');
legend(mot(1).name, mot(2).name, 'Location','northeast');

% 2) Armature current vs ω
nexttile;
plot(w, mot(1).Ia, 'LineWidth',1.8); hold on;
plot(w, mot(2).Ia, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Armature Current  I_a  [A]');
title('Armature Current vs. Speed');
legend(mot(1).name, mot(2).name, 'Location','northeast');

% 3) Armature copper loss vs ω
nexttile;
plot(w, mot(1).Ploss, 'LineWidth',1.8); hold on;
plot(w, mot(2).Ploss, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Copper Loss  I_a^2 R_a  [W]');
title('Armature Power Loss vs. Speed');
legend(mot(1).name, mot(2).name, 'Location','northeast');

% 4) Mechanical power vs ω
nexttile;
plot(w, mot(1).Pmech, 'LineWidth',1.8); hold on;
plot(w, mot(2).Pmech, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Mechanical Power  T\cdot\omega  [W]');
title('Mechanical Output Power vs. Speed');
legend(mot(1).name, mot(2).name, 'Location','northeast');

sgtitle('4.a  Ideal Motor and Battery Performance');

% Optional: save the figure
 saveas(gcf, 'hw1_q4a_ideal_motor_performance.png');

%% ================== 4.b) Non-Ideal Motor and Battery Performance ==================
% Battery non-idealities
Resr = 0.5;          % [ohm] battery ESR (series with the ideal 24V source)
Imax = 30;           % [A]   max allowable battery (armature) current

% Recompute curves with ESR and current limit
for i = 1:2
    K  = mot(i).K;
    Ra = mot(i).Ra;

    % Armature/battery current with ESR in series:
    % Vdc - Ia*Resr = K*w + Ia*Ra   =>   Ia = (Vdc - K*w) / (Ra + Resr)
    Ia = (Vdc - K.*w) ./ (Ra + Resr);

    % Enforce physical limits: no reverse current in this simple model, and cap to Imax
    Ia = max(Ia, 0);                 % no regen in this assignment's model
    Ia = min(Ia, Imax);              % battery current limit

    % Derived quantities
    T      = K .* Ia;                        % [Nm]
    PlossA = Ia.^2 .* Ra;                    % [W] armature copper loss
    PlossB = Ia.^2 .* Resr;                  % [W] battery ESR loss
    Vterm  = Vdc - Ia.*Resr;                 % [V] motor terminal voltage after ESR drop
    Pmech  = T .* w;                         % [W] mechanical output power
    Pin    = Vdc .* Ia;                      % [W] power drawn from ideal battery element
    PtoMotor = Vterm .* Ia;                  % [W] electrical power delivered to motor

    % Store for plotting/optional inspection
    mot(i).Ia_b   = Ia;
    mot(i).T_b    = T;
    mot(i).PlossA_b = PlossA;
    mot(i).PlossB_b = PlossB;
    mot(i).Pmech_b  = Pmech;
    mot(i).Pin_b    = Pin;
    mot(i).PtoMotor_b = PtoMotor;

    % Useful corner points
    mot(i).w_nl     = Vdc / K;               % [rad/s] no-load speed (unchanged by ESR if Ia=0)
    mot(i).I_stallB = min(Vdc/(Ra+Resr), Imax);
    mot(i).T_stallB = K * mot(i).I_stallB;
end

% ---------- Plots (one page) ----------
fig_b = figure('Name','4.b) Non-Ideal Motor and Battery Performance','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% (1) Torque vs speed
nexttile;
plot(w, mot(1).T_b, 'LineWidth',1.8); hold on;
plot(w, mot(2).T_b, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Torque  T  [N·m]');
title('Torque vs. Speed (with ESR & I_{max})');
legend(mot(1).name, mot(2).name, 'Location','northeast');

% (2) Armature current vs speed
nexttile;
plot(w, mot(1).Ia_b, 'LineWidth',1.8); hold on;
plot(w, mot(2).Ia_b, 'LineWidth',1.8);
yline(Imax,'k--','I_{max}','LabelVerticalAlignment','bottom');
grid on; xlabel('\omega  [rad/s]'); ylabel('Armature/Batt Current  I_a  [A]');
title('Current vs. Speed (with ESR & I_{max})');
legend(mot(1).name, mot(2).name, 'Location','northeast');

% (3) Armature copper loss vs speed
nexttile;
plot(w, mot(1).PlossA_b, 'LineWidth',1.8); hold on;
plot(w, mot(2).PlossA_b, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Armature Copper Loss  I_a^2 R_a  [W]');
title('Armature Power Loss vs. Speed (non-ideal battery)');
legend(mot(1).name, mot(2).name, 'Location','northeast');

% (4) Mechanical power vs speed
nexttile;
plot(w, mot(1).Pmech_b, 'LineWidth',1.8); hold on;
plot(w, mot(2).Pmech_b, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Mechanical Power  T\cdot\omega  [W]');
title('Mechanical Output Power vs. Speed (non-ideal battery)');
legend(mot(1).name, mot(2).name, 'Location','northeast');

sgtitle('4.b) Non-Ideal Motor and Battery Performance');

% --------- Save Figure ---------
saveas(fig_b, 'hw1_q4b_nonideal_motor_performance.png');  % saves in current folder
fprintf('Saved figure as hw1_q4b_nonideal_motor_performance.png\n');

% --------- Optional quick summary in console ---------
fprintf('\n=== 4.b) Impact of non-ideal battery (ESR=%.2f Ω, Imax=%g A) ===\n', Resr, Imax);
for i = 1:2
    fprintf('%s:\n', mot(i).name);
    fprintf('  Stall current (ideal):     %.2f A\n',  Vdc/mot(i).Ra);
    fprintf('  Stall current (with ESR):  %.2f A (capped to Imax => %.2f A)\n', Vdc/(mot(i).Ra+Resr), mot(i).I_stallB);
    fprintf('  Stall torque (with ESR/Imax): %.3f N·m\n', mot(i).T_stallB);
    fprintf('  No-load speed:             %.0f rad/s (unchanged; Ia≈0 => no ESR drop)\n', mot(i).w_nl);
end
fprintf(['Key effects: ESR reduces available voltage at high current ⇒ lower stall current and a ' ...
         'shallower torque–speed slope; current limit clamps torque at low speed; ' ...
         'no-load speed is essentially unchanged (Ia≈0 → minimal ESR drop).\n\n']);

%% ================== 4.d) Selected Motor with and without Field Weakening ==================
% Choose which motor from earlier sections to analyze:
sel = 2;   % 1 => Motor 1 (C23_L55_20), 2 => Motor 2 (C34_L60_10)

% Grab params and define weakened K
Ra   = mot(sel).Ra;
Knom = mot(sel).K;
Kweak = 0.5 * Knom;         % 50% field weakening
nameSel = mot(sel).name;

% Speed axis should cover up to the higher no-load speed (with weakened field)
w0_nom  = Vdc / Knom;
w0_weak = Vdc / Kweak;
wmax_d  = 1.05 * w0_weak;            % a little margin beyond the highest no-load speed
w_d     = linspace(0, wmax_d, 700).';

% Helper: compute curves for a given K with ESR + Imax limits
compute_curves = @(K) struct( ...
    'Ia',    min( max( (Vdc - K.*w_d) ./ (Ra + Resr), 0 ), Imax ), ...
    'K',     K );

C_nom  = compute_curves(Knom);
C_weak = compute_curves(Kweak);

% Derived quantities
T_nom   = C_nom.K  .* C_nom.Ia;
T_weak  = C_weak.K .* C_weak.Ia;

Ploss_nomA  = C_nom.Ia.^2  .* Ra;          % armature copper loss
Ploss_weakA = C_weak.Ia.^2 .* Ra;

Pmech_nom   = T_nom  .* w_d;
Pmech_weak  = T_weak .* w_d;

% ---------------- Plots ----------------
fig_d = figure('Name','4.d) Selected Motor with and without Field Weakening','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) Torque
nexttile;
plot(w_d, T_nom, 'LineWidth',1.8); hold on;
plot(w_d, T_weak, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Torque  T  [N·m]');
title(sprintf('%s — Torque vs. Speed', nameSel));
legend('No field weakening','With field weakening (K=0.5K)','Location','northeast');

% 2) Current
nexttile;
plot(w_d, C_nom.Ia, 'LineWidth',1.8); hold on;
plot(w_d, C_weak.Ia, 'LineWidth',1.8);
yline(Imax,'k--','I_{max}','LabelVerticalAlignment','bottom');
grid on; xlabel('\omega  [rad/s]'); ylabel('Armature/Batt Current  I_a  [A]');
title('Current vs. Speed (with ESR & I_{max})');
legend('No field weakening','With field weakening','Location','northeast');

% 3) Armature copper loss
nexttile;
plot(w_d, Ploss_nomA, 'LineWidth',1.8); hold on;
plot(w_d, Ploss_weakA, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Armature Copper Loss  I_a^2R_a  [W]');
title('Armature Power Loss vs. Speed');
legend('No field weakening','With field weakening','Location','northeast');

% 4) Mechanical power
nexttile;
plot(w_d, Pmech_nom, 'LineWidth',1.8); hold on;
plot(w_d, Pmech_weak, 'LineWidth',1.8);
grid on; xlabel('\omega  [rad/s]'); ylabel('Mechanical Power  T\cdot\omega  [W]');
title('Mechanical Output Power vs. Speed');
legend('No field weakening','With field weakening','Location','northeast');

sgtitle('4.d) Selected Motor with and without Field Weakening');

% -------- Save the figure --------
saveas(fig_d, 'hw1_q4d_field_weakening.png');
fprintf('Saved figure as hw1_q4d_field_weakening.png\n');

% ---------------- Max speed for a 0.5 N·m load (robust) ----------------
Tload = 0.5;  % [N·m]

% Required armature current for the requested torque
Ireq_nom  = Tload / Knom;
Ireq_weak = Tload / Kweak;

% Voltage balance with ESR: Vdc - Ia*Resr = K*w + Ia*Ra  => w = (Vdc - Ia*(Ra+Resr))/K
w_nom_raw  = (Vdc - Ireq_nom  * (Ra + Resr)) / Knom;
w_weak_raw = (Vdc - Ireq_weak * (Ra + Resr)) / Kweak;

% Feasibility checks (current limit + positive speed)
if Ireq_nom > Imax || w_nom_raw <= 0
    w_max_nom = NaN;
else
    w_max_nom = w_nom_raw;
end

if Ireq_weak > Imax || w_weak_raw <= 0
    w_max_weak = NaN;
else
    w_max_weak = w_weak_raw;
end

fprintf('\n=== 4.d) Maximum speed for %.2f N·m load (%s) ===\n', Tload, nameSel);
if isnan(w_max_nom)
    rel = Ireq_nom > Imax;  % true if current limit is the blocker
    if rel
        fprintf(' No field weakening:  infeasible (I_req=%.2f A > Imax=%g A)\n', Ireq_nom, Imax);
    else
        fprintf(' No field weakening:  infeasible (computed speed <= 0)\n');
    end
else
    fprintf(' No field weakening:  w_max = %.2f rad/s  (I_req=%.2f A)\n', w_max_nom, Ireq_nom);
end

if isnan(w_max_weak)
    rel = Ireq_weak > Imax;
    if rel
        fprintf(' With field weakening: infeasible (I_req=%.2f A > Imax=%g A)\n', Ireq_weak, Imax);
    else
        fprintf(' With field weakening: infeasible (computed speed <= 0)\n');
    end
else
    fprintf(' With field weakening: w_max = %.2f rad/s  (I_req=%.2f A)\n', w_max_weak, Ireq_weak);
end
