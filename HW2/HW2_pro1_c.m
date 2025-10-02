% File: hw2_q1d_overlay_final.m
% HW2 Q1(d): Overlay parts (a,b,c) with requested colors/styles
%  - (a) red solid
%  - (b) green, dashed for T>0 (Q1/Q2), solid for T<0 (Q3/Q4)
%  - (c) blue,  dashed for T>0 (Q1/Q2), solid for T<0 (Q3/Q4)
%  - ω = ±5000 rpm markers in purple

clear; clc;

%% Parameters
Ra      = 1.43;
K_rated = 0.0991;               % [N·m/A] = [V·s/rad]
Vmax    = 24;  Vmin = -24;
Imax    = 5;   Imin = -5;
w_rpm   = 5000;                 % rpm
w_lim   = w_rpm*2*pi/60;        % rad/s

% Colors
colA = [0.85 0.10 0.10];  % red
colB = [0.10 0.70 0.10];  % green
colC = [0.10 0.30 0.90];  % blue
colP = [0.50 0.00 0.50];  % purple

lwA = 1.8; lwB = 2.2; lwC = 2.2;

%% ---------- (a) No regen: I in [0, +Imax], K = K_rated ----------
T_a     = linspace(0, K_rated*Imax, 400);
wVmax_a =  Vmax./K_rated - (Ra./K_rated.^2).*T_a;
wVmin_a =  Vmin./K_rated - (Ra./K_rated.^2).*T_a;

%% ---------- (b) Full 4-quadrant: K = K_rated (with feasibility mask) ----------
w_b = linspace(-700, 700, 1601).';

I_low_b  = (-Vmax - K_rated*w_b)/Ra;   % from -Vmax rail
I_high_b = ( Vmax - K_rated*w_b)/Ra;   % from +Vmax rail

Imin_feas = max(I_low_b,  -Imax);      % lower feasible current
Imax_feas = min(I_high_b, +Imax);      % upper feasible current
feas = Imax_feas >= Imin_feas;         % mask where feasible

T_pos_b = K_rated .* Imax_feas;  T_pos_b(~feas) = NaN;   % +I boundary (usually T>0)
T_neg_b = K_rated .* Imin_feas;  T_neg_b(~feas) = NaN;   % –I boundary (usually T<0)

% Split by sign of torque so we can style T>0 dashed, T<0 solid
mask_b_pos = ~isnan(T_pos_b) & (T_pos_b > 0);
mask_b_neg = ~isnan(T_neg_b) & (T_neg_b < 0);

%% ---------- (c) Field-weakening to ±5000 rpm ----------
w_c = linspace(-w_lim, w_lim, 2001).';
w_eps = w_c; w_eps(abs(w_eps)<1e-9) = 1e-9;

Vrail   = @(w) Vmax*sign(w);
Kbound  = @(w,I) (Vrail(w) - I*Ra)./w;

K_up_c  = min(K_rated, Kbound(w_eps, +Imax));   % +I boundary (T_up_c ≥ 0)
K_low_c = min(K_rated, Kbound(w_eps, -Imax));   % –I boundary (T_low_c ≤ 0)

T_up_c  = (+Imax).*K_up_c;
T_low_c = (-Imax).*K_low_c;

mask_c_pos = (T_up_c > 0);   % upper half
mask_c_neg = (T_low_c < 0);  % lower half

%% ---------- Plot ----------
figure('Color','w','Name','HW2 Q1(d) — All overlays'); hold on; grid on;

% (a) RED solid (so it shows through dashed overlays)
plot(wVmax_a, T_a, 'Color',colA, 'LineWidth',lwA, 'LineStyle','-', ...
     'DisplayName','(a) K=K_{rated}, I\in[0,+I_{max}]');
plot(wVmin_a, T_a, 'Color',colA, 'LineWidth',lwA, 'LineStyle','-', ...
     'HandleVisibility','off');

% (b) GREEN — dashed in T>0 (Q1/Q2), solid in T<0 (Q3/Q4)
plot(w_b(mask_b_pos), T_pos_b(mask_b_pos), 'Color',colB, 'LineWidth',lwB, ...
     'LineStyle','--', 'DisplayName','(b) +I_{max} (T>0)');
plot(w_b(~mask_b_pos), T_pos_b(~mask_b_pos), 'Color',colB, 'LineWidth',lwB, ...
     'LineStyle','-',  'HandleVisibility','off'); % (rarely needed)

plot(w_b(mask_b_neg), T_neg_b(mask_b_neg), 'Color',colB, 'LineWidth',lwB, ...
     'LineStyle','-',  'DisplayName','(b) -I_{max} (T<0)');
plot(w_b(~mask_b_neg), T_neg_b(~mask_b_neg), 'Color',colB, 'LineWidth',lwB, ...
     'LineStyle','--', 'HandleVisibility','off'); % (rarely needed)

% (c) BLUE — dashed in T>0 (Q1/Q2), solid in T<0 (Q3/Q4)
plot(w_c(mask_c_pos), T_up_c(mask_c_pos), 'Color',colC, 'LineWidth',lwC, ...
     'LineStyle','--', 'DisplayName','(c) +I_{max} (T>0)');
plot(w_c(mask_c_neg), T_low_c(mask_c_neg), 'Color',colC, 'LineWidth',lwC, ...
     'LineStyle','-',  'DisplayName','(c) -I_{max} (T<0)');

% ω-markers in PURPLE
xline(+w_lim,'--','\omega=+5000 rpm','Color',colP,'LineWidth',1.4, ...
      'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom', ...
      'DisplayName','\omega=+5000 rpm');
xline(-w_lim,'--','\omega=-5000 rpm','Color',colP,'LineWidth',1.4, ...
      'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom', ...
      'DisplayName','\omega=-5000 rpm');

% Guides
xline(0,'k:','HandleVisibility','off'); 
yline(0,'k:','HandleVisibility','off');

xlabel('\omega  [rad/s]'); ylabel('Torque  T  [N·m]');
title('Torque–Speed Envelopes — (a)=red (solid), (b)=green, (c)=blue;  T>0 dashed, T<0 solid;  \omega=\pm5000 rpm in purple');

legend('Location','northeast');
ax = gca; ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
xlim([-700 700]); ylim([-0.55 0.55]);

saveas(gcf,'hw2_q1d_overlay_final.png');
