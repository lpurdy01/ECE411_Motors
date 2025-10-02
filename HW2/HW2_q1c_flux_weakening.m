% File: hw2_q1c_field_weakening_fixed.m
% 4-quadrant torque–speed capability with field weakening to ±5000 rpm

clear; clc;
Ra      = 1.43;          % [ohm]
K_rated = 0.0991;        % [N·m/A] == [V·s/rad]
Vmax    = 24;            % [V]
Imax    = 5;             % [A]
w_rpm   = 5000;          % [rpm]
w_lim   = w_rpm*2*pi/60; % [rad/s]

% Speed grid (exclude exact 0 to avoid 0/0; nearby points clamp to K_rated)
N = 2000;
w = linspace(-w_lim, w_lim, N).';
w(abs(w) < 1e-9) = 1e-9;

% Helper: choose the correct rail for the sign of omega (motoring direction)
Vrail = @(w) Vmax*sign(w);             % +Vmax for w>0, -Vmax for w<0

% Kmax at each (w, I) under voltage limit, then clamp to K_rated
Kmax = @(w,I) min(K_rated, (Vrail(w) - I*Ra)./w);

% Upper boundary (I = +Imax), Lower boundary (I = -Imax)
K_up   = arrayfun(@(wi) Kmax(wi, +Imax), w);
K_low  = arrayfun(@(wi) Kmax(wi, -Imax), w);

T_up   = (+Imax).*K_up;                % positive torque boundary (Q1 & Q4)
T_low  = (-Imax).*K_low;               % negative torque boundary (Q2 & Q3)

% ---- Plot ----
figure('Color','w','Name','HW2 Q1(c) Field Weakening — Fixed');
hold on; grid on;
plot(w, T_up,  'LineWidth', 2);
plot(w, T_low, 'LineWidth', 2);
xline(0,'k:'); yline(0,'k:');
xlabel('\omega  [rad/s]'); ylabel('Torque  T  [N·m]');
title('Torque–Speed Capability with Field Weakening to \pm5000 rpm');
legend('+I_{max} boundary','-I_{max} boundary','Location','northeast');

% ---- Report the four required points at ±5000 rpm ----
cases = {
    'ω=+5000 rpm, +T (Q1 motoring)',  +w_lim, +Imax
    'ω=+5000 rpm, -T (Q2 generating)',+w_lim, -Imax
    'ω=-5000 rpm, +T (Q4 generating)',-w_lim, +Imax
    'ω=-5000 rpm, -T (Q3 motoring)',  -w_lim, -Imax
};

fprintf('=== HW2 Q1(c): outputs at ±%d rpm ===\n', w_rpm);
for k = 1:size(cases,1)
    label = cases{k,1};  w0 = cases{k,2};  I = cases{k,3};

    % compute K*, T, Vs at that point
    Kstar = Kmax(w0, I);
    T     = I*Kstar;
    Vs    = Kstar*w0 + I*Ra;      % should be sign(w0)*Vmax
    Pmech = T*w0;                 % signed shaft power

    % "Output power" convention: mechanical out for motoring, electrical out for regen
    if sign(T)==sign(w0)          % motoring
        P_out = Pmech;            % mechanical output (can be negative if w0<0, but magnitude is what matters)
    else                          % generating
        P_out = abs(Pmech);       % electrical output magnitude
    end

    fprintf('%s:\n', label);
    fprintf('  K = %.6f N·m/A\n', Kstar);
    fprintf('  T = %.4f N·m\n', T);
    fprintf('  P_out = %.3f kW\n', P_out/1000);
    fprintf('  V_s = %.2f V\n\n', Vs);
end
