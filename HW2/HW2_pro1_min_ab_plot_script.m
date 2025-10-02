% File: hw2_q1_capability_plots_only.m
% Purpose: Plot torque–speed capability envelopes for HW2 Q1(a) and Q1(b)
% Motor: C23_L55_20

clear; clc;

% -------- Motor params --------
Ra = 1.43;            % [ohm]
K  = 0.0991;          % [N·m/A] == [V·s/rad]
Vmax = 24; Vmin = -24;

% Helper to plot one envelope
plot_envelope = @(Imin, Imax, ttl, figname) ...
    local_plot_envelope(Imin, Imax, Vmin, Vmax, Ra, K, ttl, figname);

% ---- (a) No regen: I in [0, +5] ----
plot_envelope(0, 5, ...
    'HW2 Q1(a): Capability (I \in [0,+5] A, V \in [-24,+24] V)', ...
    'hw2_q1a_capability.png');

% ---- (b) Full 4-quadrant: I in [-5, +5] ----
plot_envelope(-5, 5, ...
    'HW2 Q1(b): Capability (I \in [-5,+5] A, V \in [-24,+24] V)', ...
    'hw2_q1b_capability.png');

% ================== helpers ==================
function local_plot_envelope(Imin, Imax, Vmin, Vmax, Ra, K, ttl, figname)
    % Boundary segments on T–ω plane are given by voltage lines:
    % ω = Va/K - (Ra/K^2) * T  with Va ∈ {Vmin, Vmax}
    % and by current limits T = K*I with I ∈ [Imin, Imax]

    % Build a torque axis that covers both current limits
    T = linspace(K*Imin, K*Imax, 400);

    % Voltage-limit lines (upper and lower)
    w_Vmax =  Vmax./K - (Ra./K.^2).*T;
    w_Vmin =  Vmin./K - (Ra./K.^2).*T;

    % Current-limit lines (horizontal)
    T_Imax =  K*Imax;
    T_Imin =  K*Imin;   % equals 0 in (a)

    % Compute their ω extents at the current limits
    w_Imax_Vmax = (Vmax - Imax*Ra)/K;
    w_Imax_Vmin = (Vmin - Imax*Ra)/K;
    w_Imin_Vmax = (Vmax - Imin*Ra)/K;
    w_Imin_Vmin = (Vmin - Imin*Ra)/K;

    % Figure
    figure('Color','w','Name',ttl); hold on; grid on;
    xlabel('\omega  [rad/s]'); ylabel('Torque  T  [N·m]'); title(ttl);

    % Voltage lines
    plot(w_Vmax, T, 'LineWidth',1.8);
    plot(w_Vmin, T, 'LineWidth',1.8);

    % Current-limit segments (horizontal)
    plot([w_Imax_Vmin, w_Imax_Vmax],[T_Imax, T_Imax],'k','LineWidth',2);
    plot([w_Imin_Vmin, w_Imin_Vmax],[T_Imin, T_Imin],'k--','LineWidth',1.2);

    % Axes crossing
    xline(0,'k:'); yline(0,'k:');

    legend('V = +V_{max}','V = -V_{max}', ...
           'I = +I_{max}','I = I_{min}','Location','best');
    saveas(gcf, figname);
end
