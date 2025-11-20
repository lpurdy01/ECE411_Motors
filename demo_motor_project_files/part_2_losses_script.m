%% clear all close all section
clear all; close all; clc;

%% AI disclaimer
% The following code was generated with the assistance of GPT-5 and 
% GPT-5-Codex from GitHub Copilot.
% It has been based on existing human written code,
% human-reviewed, and progressively refined.

%% Init section
% Script settings
outputs_directory = 'outputs';
if ~exist(outputs_directory, 'dir')
    mkdir(outputs_directory);
end

save_graphs = true; % set to true to save graphs to outputs directory

P = machine_params(); % Load machine parameters

% --- field current cases ---
field_current_cases = [7.0 4.666 2.333]; % A dc, cases i, ii, iii

% handy constants
V0_limit = P.Vdc/sqrt(3);          % linear modulation ceiling (LN peak)
I0_limit = P.I0_max;               % stator current peak limit
P_coreR  = 500;                    % W, provided reference point
k_inv    = 0.015;                  % inverter loss coefficient
E0_R     = P.wb * P.Lm * P.If_max; % EMF magnitude at rated field & base freq

fprintf('Loss model constants:\n');
fprintf('  V0 limit (peak LN) = %.1f V\n', V0_limit);
fprintf('  I0 limit (peak)    = %.1f A\n', I0_limit);
fprintf('  P_coreR             = %.1f W\n', P_coreR);
fprintf('  k_inv               = %.4f\n', k_inv);
fprintf('  E0_R                = %.1f V\n\n', E0_R);

%% Simplifying assumptions and notes:
%
% 1. Id = 0 (all-torque current) and sinusoidal steady-state.
% 2. Losses are computed directly; their secondary influence on currents
%    and voltages is ignored per instructions.
% 3. Magnetic model remains linear with constant parameters.


%% --- Plot 1: Field winding I^2R ---
If_axis = linspace(0, P.If_max, 400);
Pf_axis = (If_axis.^2) * P.Rf;

fig = figure; hold on; grid on;
plot(If_axis, Pf_axis, 'LineWidth', 2);
xlabel('Field current  I_f  (A)');
ylabel('Field loss  P_f  (W)');
title('Field Winding I^2R Loss vs. Field Current');

if save_graphs && isvalid(fig)
    saveas(fig, fullfile(outputs_directory, 'part_2_loss_plot_field.png'));
end


%% --- Plot 2: Core loss vs EMF magnitude ---
E0_axis = linspace(0, E0_R, 400);
Pcore_axis = P_coreR * (E0_axis / E0_R);

fig = figure; hold on; grid on;
plot(E0_axis, Pcore_axis, 'LineWidth', 2);
xlabel('Back-EMF magnitude  E_0  (V_{LN,peak})');
ylabel('Core loss  P_{core}  (W)');
title('Core Loss vs. Back-EMF Magnitude');

if save_graphs && isvalid(fig)
    saveas(fig, fullfile(outputs_directory, 'part_2_loss_plot_core.png'));
end


%% --- Plot 3: Stator+inverter loss vs stator current ---
I0_axis = linspace(0, P.I0_max, 400);
Ps_axis   = 0.5 * (I0_axis.^2) * P.Rs;
Pinv_axis = k_inv * P.Vdc * I0_axis;
Ps_total_axis = Ps_axis + Pinv_axis;

fig = figure; hold on; grid on;
plot(I0_axis, Ps_total_axis, 'LineWidth', 2);
xlabel('Stator current magnitude  |I_0|  (A peak)');
ylabel('Loss  (W)');
title('Stator I^2R + Inverter Loss vs. Stator Current');

if save_graphs && isvalid(fig)
    saveas(fig, fullfile(outputs_directory, 'part_2_loss_plot_stator_inverter.png'));
end




% Speed and torque sweeps
mechanical_velocity_om_m = linspace(0, 3*P.wmb, 200);  % mechanical rad/s
Te_values_Nm = linspace(0, 220, 200);                  % torque range

[Omega_m_grid_SI, Te_grid_SI] = meshgrid(mechanical_velocity_om_m, Te_values_Nm);

for k = 1:numel(field_current_cases)
    If = field_current_cases(k);
    flux_linkage = P.Lm * If; % constant for this case

    % Electrical speed grid and air-gap EMF
    Omega_e_grid_SI = (P.p/2) .* Omega_m_grid_SI;
    Eq_grid = Omega_e_grid_SI * flux_linkage;

    % Stator current needed for requested torque (Id = 0)
    Iq_grid = (4/(3*P.p)) * (Te_grid_SI / flux_linkage);
    I_mag_grid = abs(Iq_grid);

    % Voltage requirements (same KVL used previously)
    Vq_grid = Iq_grid * P.Rs + Eq_grid;
    Vd_grid = -Omega_e_grid_SI .* P.Ls .* Iq_grid;
    V0_grid = sqrt(Vq_grid.^2 + Vd_grid.^2);

    current_mask = I_mag_grid <= I0_limit;
    voltage_mask = V0_grid <= V0_limit;
    feasible_mask = current_mask & voltage_mask;

    % Loss components
    Pf_scalar = If^2 * P.Rf;
    E0_grid = abs(Eq_grid);
    Pcore_grid = P_coreR * (E0_grid / E0_R);
    Ps_grid = 0.5 * (I_mag_grid.^2) * P.Rs;
    Pinv_grid = k_inv * P.Vdc * I_mag_grid;

    P_total_grid = Pf_scalar + Pcore_grid + Ps_grid + Pinv_grid;
    P_total_grid(~feasible_mask) = NaN;

    % ---- Plot
    fig = figure; hold on; grid on;
    contourf(Omega_m_grid_SI, Te_grid_SI, P_total_grid/1e3, 20, 'ShowText', true);
    cb = colorbar;
    ylabel(cb, 'Total loss  P_{loss}  (kW)');
    xlabel('\omega_m  (mechanical rad/s)');
    ylabel('Electromagnetic torque  T_e  (N\cdotm)');
    title(sprintf('Total Power Loss Contours  |  I_f = %.3f A', If));
    xlim([mechanical_velocity_om_m(1), mechanical_velocity_om_m(end)]);
    ylim([Te_values_Nm(1), Te_values_Nm(end)]);
    caxis([0, 5]);

    % Save graph if enabled
    if save_graphs && isvalid(fig)
        filename = sprintf('part_2_loss_contours_If_%0.3fA.png', If);
        saveas(fig, fullfile(outputs_directory, filename));
    end
end


%% Operating scenarios and recommendations
scenarios = struct( ...
    'name', {'Near max acceleration', 'Stop-and-go city traffic', 'Highway cruising'}, ...
    'omega_m', {250, 500, 1250}, ...
    'torque', {200, 60, 15});

num_scenarios = numel(scenarios);
num_components = 4; % Pf, Pcore, Ps, Pinv
loss_components = nan(num_scenarios, num_components);
selected_field_currents = nan(num_scenarios, 1);

fprintf('Scenario analysis (totals in watts):\n');
for s = 1:num_scenarios
    scenario = scenarios(s);
    best_total = inf;
    best_idx = NaN;
    best_components = nan(1, num_components);

    for k = 1:numel(field_current_cases)
        If = field_current_cases(k);
        flux_linkage = P.Lm * If;

        if abs(flux_linkage) < 1e-9
            continue; % no excitation, cannot sustain torque
        end

        omega_e = (P.p/2) * scenario.omega_m;
        Eq = omega_e * flux_linkage;

        Iq = (4/(3*P.p)) * (scenario.torque / flux_linkage);
        I_mag = abs(Iq);
        Vq = Iq * P.Rs + Eq;
        Vd = -omega_e * P.Ls * Iq;
        V0 = sqrt(Vq^2 + Vd^2);

        if I_mag > I0_limit || V0 > V0_limit
            continue; % infeasible operating point
        end

    Pf = If^2 * P.Rf;
        E0 = abs(Eq);
        Pcore = P_coreR * (E0 / E0_R);
        Ps = 0.5 * (I_mag^2) * P.Rs;
        Pinv = k_inv * P.Vdc * I_mag;
        total_loss = Pf + Pcore + Ps + Pinv;

        if total_loss < best_total
            best_total = total_loss;
            best_idx = k;
            best_components = [Pf, Pcore, Ps, Pinv];
        end
    end

    if isnan(best_idx)
        fprintf('  %-24s : no feasible field current found.\n', scenario.name);
    else
        selected_field_currents(s) = field_current_cases(best_idx);
        loss_components(s, :) = best_components;
        fprintf('  %-24s : If = %.3f A,   P_loss = %.1f W\n', ...
            scenario.name, selected_field_currents(s), best_total);
    end
end


%% Extra credit: loss breakdown bar chart
valid_rows = all(isfinite(loss_components), 2);
if any(valid_rows)
    scenario_names = {scenarios(valid_rows).name};
    fig = figure; hold on; grid on;
    bar(loss_components(valid_rows, :) / 1e3, 'stacked');
    set(gca, 'XTickLabel', scenario_names, 'XTickLabelRotation', 15);
    ylabel('Loss contribution (kW)');
    title('Loss Breakdown by Scenario (recommended I_f per case)');
    legend({'P_f', 'P_{core}', 'P_s', 'P_{inv}'}, 'Location', 'northwest');

    if save_graphs && isvalid(fig)
        saveas(fig, fullfile(outputs_directory, 'part_2_loss_breakdown_bar.png'));
    end
end

fprintf('\nSection 2 loss analysis complete.\n');
