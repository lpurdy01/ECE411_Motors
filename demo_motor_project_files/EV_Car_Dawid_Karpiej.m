clear; close all;

% machine parameters

Lm  = 17.5e-3;
Ls  = 0.6e-3;
Rs  = 0.0155;
Rf  = 26.5;
p   = 8;
wb  = 2*pi*267;
wmb = 2*wb/p;

% inverter

Vdc    = 550;
I0_max = 300;

% field current cases

If_i = 7.0;
If_ii = 4.666;
If_iii = 2.333;

% flux cases

flux_i = Lm * If_i;
flux_ii = Lm * If_ii;
flux_iii = Lm * If_iii;

% PART 1 -----------------------------------------------------------------

% PLOT 1 -----------------------------------------------------------------

% data points for x axis

wm_b = linspace(0, 3*wmb, 1000);

% Separate q and d components

Eq_i = (wm_b*p/2) * Lm * If_i;
Eq_ii = (wm_b*p/2) * Lm * If_ii;
Eq_iii = (wm_b*p/2) * Lm * If_iii;
Ed = 0;

Iq = I0_max;
Id = 0;

% general KVL solution for SM with nonzero resistance

Vq_i = Iq*Rs + (wm_b*p/2)*Ls*Id + Eq_i;
Vd_i = Id*Rs - (wm_b*p/2)*Ls*Iq + Ed;

Vq_ii = Iq*Rs + (wm_b*p/2)*Ls*Id + Eq_ii;
Vd_ii = Id*Rs - (wm_b*p/2)*Ls*Iq + Ed;

Vq_iii = Iq*Rs + (wm_b*p/2)*Ls*Id + Eq_iii;
Vd_iii = Id*Rs - (wm_b*p/2)*Ls*Iq + Ed;

% voltage magnitude

V0_i = sqrt(Vq_i.^2 + Vd_i.^2);
V0_ii = sqrt(Vq_ii.^2 + Vd_ii.^2);
V0_iii = sqrt(Vq_iii.^2 + Vd_iii.^2);

% plot

figure
plot(wm_b, V0_i, wm_b, V0_ii, wm_b, V0_iii)
title('Plot 1 Voltage Amplitude vs Mechanical Speed')
xlabel('Speed [rad/s]')
ylabel('Voltage Amplitude')
legend('If = 7.0', 'If = 4.6', 'If = 2.3')

% Pick out the points for each If case where V0 = Vdc/sqrt(3)
V0_limit = Vdc/sqrt(3);

[~, idx_i] = min(abs(V0_i - V0_limit));
cross_speed_i = wm_b(idx_i);

[~, idx_ii] = min(abs(V0_ii - V0_limit));
cross_speed_ii = wm_b(idx_ii);

[~, idx_iii] = min(abs(V0_iii - V0_limit));
cross_speed_iii = wm_b(idx_iii);

% print the Outputs
fprintf('Crossing Speeds for V0 = Vdc/sqrt(3):\n')
fprintf('If = 7.0 A :  omega_m* = %.1f rad/s  (%.0f rpm)\n', ...
    cross_speed_i, cross_speed_i*60/(2*pi));
fprintf('If = 4.6 A :  omega_m* = %.1f rad/s  (%.0f rpm)\n', ...
    cross_speed_ii, cross_speed_ii*60/(2*pi)); 
fprintf('If = 2.3 A :  omega_m* = %.1f rad/s  (%.0f rpm)\n', ...
    cross_speed_iii, cross_speed_iii*60/(2*pi));



% PLOT 2 -----------------------------------------------------------------

% data points for x and y axes

Iq = linspace(0, 300, 1000);
wm_b = linspace(0, 3*wmb, 1000);

Iq_col = Iq.';               % column vector for outer products
we_row = (p/2) * wm_b;        % row vector of electrical speed
Vd_grid_common = -(Iq_col * we_row) * Ls;
Vq_base = Iq_col * Rs;

% torque equation

Te_i = 3*p/4*flux_i*Iq;
Te_ii = 3*p/4*flux_ii*Iq;
Te_iii = 3*p/4*flux_iii*Iq;

% mechanical power

Pm_i = Te_i'*wm_b;
Pm_ii = Te_ii'*wm_b;
Pm_iii = Te_iii'*wm_b;

Vq_i_grid = Vq_base + Eq_i; % copper losses + bemf
V0_i_grid = sqrt(Vq_i_grid.^2 + Vd_grid_common.^2); % magnitude
Pm_i(V0_i_grid > V0_limit) = NaN; % check if unfeasible

Vq_ii_grid = Vq_base + Eq_ii;
V0_ii_grid = sqrt(Vq_ii_grid.^2 + Vd_grid_common.^2);
Pm_ii(V0_ii_grid > V0_limit) = NaN;

Vq_iii_grid = Vq_base + Eq_iii;
V0_iii_grid = sqrt(Vq_iii_grid.^2 + Vd_grid_common.^2);
Pm_iii(V0_iii_grid > V0_limit) = NaN;

% plot

figure
contourf(wm_b, Te_i, Pm_i, 10, "ShowText", true)
title('Plot 2 Mechanical Power vs Mech Speed and Torque with I_f=7.0')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Mechanical Power')
ylim([0 220])
colorbar

figure
contourf(wm_b, Te_ii, Pm_ii, 10, "ShowText", true)
title('Plot 2 Mechanical Power vs Mech Speed and Torque with I_f=4.6')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Mechanical Power')
ylim([0 220])
colorbar

figure
contourf(wm_b, Te_iii, Pm_iii, 10, "ShowText", true)
title('Plot 2 Mechanical Power vs Mech Speed and Torque with I_f=2.3')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Mechanical Power')
ylim([0 220])
colorbar

% PLOT 3 -----------------------------------------------------------------

% data points for x and y

wm_b = linspace(0, 3*wmb, 1000);
Iq = linspace(0, 300, 1000);
Id = 0;

% calculate vq and vd for varying current and speed

Vq_i = Iq*Rs + Eq_i;
Vd_i = (wm_b*p/2)'*Ls.*Iq;

Vq_ii = Iq*Rs + Eq_ii;
Vd_ii = (wm_b*p/2)'*Ls.*Iq;

Vq_iii = Iq*Rs + Eq_iii;
Vd_iii = (wm_b*p/2)'*Ls.*Iq;

% calculate modulation in terms of vd and vq

Mq_i = 2*Vq_i/Vdc;
Mq_ii = 2*Vq_ii/Vdc;
Mq_iii = 2*Vq_iii/Vdc;

Md_i = 2*Vd_i/Vdc;
Md_ii = 2*Vd_ii/Vdc;
Md_iii = 2*Vd_iii/Vdc;

% modulation amplitude

M0_i = sqrt(Mq_i.^2 + Md_i.^2);
M0_ii = sqrt(Mq_ii.^2 + Md_ii.^2);
M0_iii = sqrt(Mq_iii.^2 + Md_iii.^2);

% plot 

figure
contourf(wm_b, Te_i, M0_i, 10, "ShowText", true)
title('Plot 3 Modulatiuon Amplitude vs Mechanical Speed and Torque with I_f=7.0')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Modulation Index')
colorbar

figure
contourf(wm_b, Te_ii, M0_ii, 10, "ShowText", true)
title('Plot 3 Modulatiuon Amplitude vs Mechanical Speed and Torque with I_f=4.6')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Modulation Index')
colorbar

figure
contourf(wm_b, Te_iii, M0_iii, 10, "ShowText", true)
title('Plot 3 Modulatiuon Amplitude vs Mechanical Speed and Torque with I_f=2.3')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Modulation Index')
colorbar

% PART 2 -----------------------------------------------------------------

% PLOT 1 -----------------------------------------------------------------

% calculate field losses for all cases, constant for any speed or torque

If = linspace(0, 7, 1000);

Pf = (If.^2)*Rf;

% for total power calculation, use all three cases

Pf_i = (If_i^2)*Rf;
Pf_ii = (If_ii^2)*Rf;
Pf_iii = (If_iii^2)*Rf;

% plot

figure
plot(If, Pf)
title('Plot 1 Field Winding Loss vs Field Current')
xlabel('Field Current')
ylabel('Field Winding Loss')

% PLOT 2 -----------------------------------------------------------------

% calculate magnetic core losses, varying with speed

P_r = 500;
E0_r = (wb)*Lm*If_i;

E0 = linspace(0, E0_r, 1000);

Pcore = P_r*E0/E0_r;

% for total power calculation, use all three cases

Pcore_i = P_r*Eq_i/E0_r;
Pcore_ii = P_r*Eq_ii/E0_r;
Pcore_iii = P_r*Eq_iii/E0_r;

% plot

figure
plot(E0, Pcore)
title('Plot 2 Core Loss vs Back EMF Magnitude')
xlabel('Back EMF Magnitude')
ylabel('Core Loss')

% PLOT 3 -----------------------------------------------------------------

% calculate stator winding losses, varying with torque

Ps = 1/2*(abs(Iq).^2)*Rs;

% inverter losses

k_inv = 0.015;

Pinv = k_inv*Vdc*Iq;

% plot

figure
plot(Iq, Ps + Pinv)
title('Plot 3 Stator and Inverter Losses vs Stator Current')
xlabel('Current Magnitude')
ylabel('Stator and Inverter Loss')

% PLOT 4 -----------------------------------------------------------------

% prepare matrices for addition by replicating matrices

Pcore_i = repmat(Pcore_i, 1000, 1);
Pcore_ii = repmat(Pcore_ii, 1000, 1);
Pcore_iii = repmat(Pcore_iii, 1000, 1);

% replicate and transpose the matries that vary with speed

Ps = transpose(repmat(Ps, 1000, 1));

Pinv = transpose(repmat(Pinv, 1000, 1));

% total power loss

Ptot_i = Pcore_i + Ps + Pinv + Pf_i;
Ptot_ii = Pcore_ii + Ps + Pinv + Pf_ii;
Ptot_iii = Pcore_iii + Ps + Pinv + Pf_iii;

Ptot_i(V0_i_grid > V0_limit) = NaN;
Ptot_ii(V0_ii_grid > V0_limit) = NaN;
Ptot_iii(V0_iii_grid > V0_limit) = NaN;

% plot

figure
contourf(wm_b, Te_i, Ptot_i, 10, "ShowText", true)
ylim([0 220])
xlim([0 3*wmb])
caxis([0 5000])
title('Plot 4 Total Power Loss vs Mechanical Speed and Torque with I_f=7.0')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Total Power Loss')
colorbar

figure
contourf(wm_b, Te_ii, Ptot_ii, 10, "ShowText", true)
ylim([0 220])
xlim([0 3*wmb])
caxis([0 5000])
title('Plot 4 Total Power Loss vs Mechanical Speed and Torque with I_f=4.6')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Total Power Loss')
colorbar

figure
contourf(wm_b, Te_iii, Ptot_iii, 10, "ShowText", true)
ylim([0 220])
xlim([0 3*wmb])
caxis([0 5000])
title('Plot 4 Total Power Loss vs Mechanical Speed and Torque with I_f=2.3')
xlabel('Mechanical Speed [rad/s]')
ylabel('Electromagnetic Torque')
zlabel('Total Power Loss')
colorbar
