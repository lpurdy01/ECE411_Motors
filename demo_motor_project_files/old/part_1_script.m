%% clear all close all section
clear all; close all; clc;

% Script settings
outputs_directory = 'outputs';
if ~exist(outputs_directory, 'dir')
    mkdir(outputs_directory);
end

save_graphs = true; % set to true to save graphs to outputs directory

%% --- machine paramters ---
Lm  = 17.5e-3;     % magnetizig induct H
Ls  = 0.6e-3;    % stator induct (q=d) [H]
Rs  = 0.0155;      % stator resist ohm
Rf  = 26.5;      % field resist [ohm]
p   = 8;          % poles
wb  = 2*pi*267;    % elec base w rad/s
wmb = 419;        % mech base w rad/s

% --- inverter / field ctrl ---
Vdc    = 600;      % dc bus V
I0_max = 300;   % max ac curr (phase peak) A
If_max = 7.0;      % max field curr A dc

% --- derived ---
Xs_at_elec_base = wb*Ls;                  % sync react ohm
Xm_at_elec_base = wb*Lm;                % mag react ohm
f_e_at_elect_base = wb/(2*pi);             % elec base freq hz
n_sync_rpm_at_elect_base = wmb*60/(2*pi);  % mech base rpm

% --- field current cases --
field_current_cases = [7.0 4.666 2.333]; % A dc, cases i, ii, iii

%% --- print summary ---
fprintf('machine params loaded:\n')
fprintf('  Lm   = %.4f H\n', Lm)
fprintf('  Ls   = %.4f H\n', Ls)
fprintf('  Rs   = %.5f ohm\n', Rs)
fprintf('  Rf   = %.2f ohm\n', Rf)
fprintf('  p    = %d poles\n', p)
fprintf('  wb   = %.1f rad/s  (%.1f hz)\n', wb, f_e_at_elect_base)
fprintf('  wmb  = %.1f rad/s  (%.1f rpm)\n', wmb, n_sync_rpm_at_elect_base)
fprintf('\ninverter params:\n')
fprintf('  Vdc     = %.0f V\n', Vdc)
fprintf('  I0_max  = %.0f A pk\n', I0_max)
fprintf('  If_max  = %.1f A dc\n', If_max)

%% Simplifying assumptions and notes:

% 1. All inverter current is used for Torque (no AC field weakening), solve KVL from there… 
% 2. Machine parameters are linear and constant (no saturation, no heating of resistances)
% 3. Don't bother with Svpwm or complex offset modulation, just go straight linear modulation 


%% --- Plot 1 ------
%
%██████  ██       ██████  ████████      ██ 
%██   ██ ██      ██    ██    ██        ███ 
%██████  ██      ██    ██    ██         ██ 
%██      ██      ██    ██    ██         ██ 
%██      ███████  ██████     ██         ██ 
%                                          
                                          
% Use https://patorjk.com/software/taag/#p=display&f=ANSI+Regular&t=PLOT+1&x=none&v=4&h=4&w=80&we=false
% for ASCII art title

% Speed sweep (mechanical rad/s); cover up to 2*wmb per spec
mechanical_velocity_om_m = linspace(0, 2*wmb, 500);  % mechanical rad/s
electrical_angular_velocity = (p/2).*mechanical_velocity_om_m;             % electrical rad/s

% Inverter linear-modulation limit
% This comes from M_0 = V0/(0.5*Vdc) <= 2/sqrt(3)  therefore  V0 <= Vdc/sqrt(3)
V0_lim = Vdc/sqrt(3);

% KVL around SM Generic Case:  V_qd = E_qd + (Rs + j*ωe*Ls) * I_qd
% or V_qd = E_qd + Zs * I_qd
% or:  I_qd = (V_qd - E_qd) / (Rs + j*ωe*Ls)) (L2-8 slide 3)
% We fix |I| = I0_max with all-torque current: I = I0_max + j*0 (Id=0) for max torque no phase.

V0_all = zeros(numel(field_current_cases), numel(mechanical_velocity_om_m));
Vq_all = zeros(size(V0_all));
Vd_all = zeros(size(V0_all));
cross_idx   = nan(size(field_current_cases));
cross_speed = nan(size(field_current_cases));

for k = 1:numel(field_current_cases)
    If = field_current_cases(k); % grab our field current case

    % Electrical speed array, impedance and EMF phasors
    Zs = Rs + 1i*electrical_angular_velocity*Ls;  % stator impedance vs speed (L2-8 slide 4)
    E  = (electrical_angular_velocity * Lm * If) + 1i*0;  % E_qd phasor: Eq=ωe*Lm*If, Ed=0
    % Using an E aligned phasor version

    % Current phasor at limit (all torque current => Id=0)
    I = I0_max + 1i*0;  % I_qd phasor: Iq=I0_max, Id=0

    % Required stator voltage phasor (from KVL)
    V = E + Zs.*I;  % V_qd phasor vs speed

    % Break into LN peak q–d components (Lipo–Novotny negative imag d-axis)
    Vq = real(V);
    Vd = -imag(V);

    % Phase fundamental (LN peak) amplitude
    V0 = abs(V);

    % Store
    Vq_all(k,:) = Vq;
    Vd_all(k,:) = Vd;
    V0_all(k,:) = V0;

    % First crossing of linear-modulation limit
    idx = find(V0 >= V0_lim, 1, 'first');
    if ~isempty(idx)
        cross_idx(k)   = idx;
        cross_speed(k) = mechanical_velocity_om_m(idx);
    end
end

%% Plotting Section
% This plotting section does use AI to get nice symbols and formatting to make plots better.
% Use ChatGPT GPT-5 and gh:copilot GPT-5-Codex

% ---- Plot
figure; hold on; grid on;
colors = lines(numel(field_current_cases));
for k = 1:numel(field_current_cases)
    plot(mechanical_velocity_om_m, V0_all(k,:), 'LineWidth', 2, 'Color', colors(k,:));
    if ~isnan(cross_idx(k))
        plot(cross_speed(k), V0_lim, 'o', 'Color', colors(k,:), ...
            'MarkerFaceColor', colors(k,:), 'MarkerSize', 7);
    end
end
yline(V0_lim, '--k', 'LineWidth', 1.5, 'Label', sprintf('V_0 limit = V_{dc}/\\surd3 = %.1f V', V0_lim), ...
    'LabelVerticalAlignment', 'bottom');

xlabel('\omega_m  (mechanical rad/s)');
ylabel('Phase fundamental amplitude  V_0  (V_{LN,peak})');
title('Inverter V_0 vs. Mechanical Speed at I_0 = I_{0,max}');
legend(arrayfun(@(x)sprintf('I_f = %.1f A',x), field_current_cases, 'UniformOutput',false), ...
       'Location', 'northwest');

% Report crossing speeds
fprintf('\n--- V0 limit crossings (V0 = Vdc/sqrt(3)) ---\n');
for k = 1:numel(field_current_cases)
    if ~isnan(cross_idx(k))
        fprintf('If = %.1f A :  omega_m* = %.1f rad/s  (%.0f rpm)\n', ...
            field_current_cases(k), cross_speed(k), cross_speed(k)*60/(2*pi));
    else
        fprintf('If = %.1f A :  no crossing in sweep range.\n', field_current_cases(k));
    end
end

% Save graph if save_graphs is true
if save_graphs
    saveas(gcf, fullfile(outputs_directory, 'part_1_plot_1_inverter_V0_vs_speed.png'));
end

