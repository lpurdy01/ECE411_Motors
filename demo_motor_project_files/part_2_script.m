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

% --- choose per-unit bases (match Part 1 assumptions) ---
VB_LN_peak = P.Vdc/sqrt(3);               % inverter linear-modulation ceiling (LN peak)
VB_LN_rms  = VB_LN_peak/sqrt(2);
PB_3ph     = 1.5 * VB_LN_peak * P.I0_max;  % ensures I_base equals I0_max (peak)
B          = pu_bases(VB_LN_rms, PB_3ph, P.wb, P.p);

fprintf('per-unit bases (LN peak):\n');
fprintf('  V_B = %.1f V\n', B.VB);
fprintf('  I_B = %.1f A\n', B.IB);
fprintf('  Z_B = %.3f ohm\n', B.ZB);
fprintf('  P_B = %.1f kW\n', B.PB/1e3);
fprintf('  T_B = %.1f N*m\n\n', B.TB);

%% Simplifying assumptions and notes:

% 1. All inverter current is used for Torque (no AC field weakening), solve KVL from there…
% 2. Machine parameters are linear and constant (no saturation, no heating of resistances)
% 3. Don't bother with Svpwm or complex offset modulation, just go straight linear modulation


%% --- Plot 2 ------
%
%██████  ██       ██████  ████████     ██████  
%██   ██ ██      ██    ██    ██             ██ 
%██████  ██      ██    ██    ██         █████  
%██      ██      ██    ██    ██        ██      
%██      ███████  ██████     ██        ███████ 
%                                          


% Speed sweep (mechanical rad/s); cover up to 2*wmb per spec
mechanical_velocity_om_m = linspace(0, 3*P.wmb, 200);  % mechanical rad/s
mechanical_velocity_pu   = mechanical_velocity_om_m / B.wmB;

% Torque sweep (N*m); spec says 0 to 230 for all plots
Te_values_Nm = linspace(0, 230, 200);
Te_values_pu = Te_values_Nm / B.TB;

% Meshgrids for contour evaluation
[Omega_m_grid_pu, Te_grid_pu] = meshgrid(mechanical_velocity_pu, Te_values_pu);
[Omega_m_grid_SI, Te_grid_SI] = meshgrid(mechanical_velocity_om_m, Te_values_Nm);

% Limits in per-unit
I_lim_pu   = P.I0_max / B.IB;
V0_lim     = P.Vdc/sqrt(3);
V0_lim_pu  = V0_lim / B.VB;

% Preallocate storage for reporting (optional)
max_torque_feasible = zeros(size(field_current_cases));

for k = 1:numel(field_current_cases)
    If = field_current_cases(k);

    % Electrical speed (SI and per-unit) across the grid
    Omega_e_grid_SI = (P.p/2) .* Omega_m_grid_SI;
    Omega_e_grid_pu = Omega_e_grid_SI / B.wB;

    % Internal EMF in per-unit (E_q component)
    Eq_grid_pu = (Omega_e_grid_SI .* P.Lm .* If) / B.VB;

    % Avoid divide-by-zero at zero speed / zero excitation by marking as undefined
    small_emf_mask = abs(Eq_grid_pu) < 1e-9;
    Eq_grid_pu(small_emf_mask) = NaN;

    % Compute required Iq in per-unit from torque definition (Id=0)
    Iq_grid_pu = (Te_grid_pu .* Omega_e_grid_pu) ./ Eq_grid_pu;
    Iq_grid_pu(~isfinite(Iq_grid_pu)) = NaN;
    Id_grid_pu = zeros(size(Iq_grid_pu));

    % Current magnitude (since Id = 0, |I| = |Iq|)
    I_mag_pu = abs(Iq_grid_pu);
    current_mask = I_mag_pu <= I_lim_pu;

    % Convert currents back to SI for helper function (expects A peak)
    Iq_grid_SI = Iq_grid_pu * B.IB;
    Id_grid_SI = Id_grid_pu * B.IB;

    % Flatten for function call and reshape results
    Iq_vec = Iq_grid_SI(:);
    Id_vec = Id_grid_SI(:);
    Om_vec = Omega_m_grid_SI(:);

    [~, Vq_pu_vec, Vd_pu_vec, V0_pu_vec, ~, ~, Eq_pu_vec] = sm_required_voltage_park_pu( ...
        Om_vec, If, Iq_vec, Id_vec, P, B);

    V0_pu_grid = reshape(V0_pu_vec, size(Omega_m_grid_pu));
    Eq_pu_grid = reshape(Eq_pu_vec, size(Omega_m_grid_pu));

    voltage_mask = V0_pu_grid <= V0_lim_pu;

    % Mechanical power in pu and SI (Id=0 assumption)
    Pmech_pu = real(Eq_pu_grid) .* Iq_grid_pu;
    Pmech_SI = Pmech_pu * B.PB; % watts

    % Apply feasibility masks
    feasible_mask = current_mask & voltage_mask;
    Pmech_SI(~feasible_mask) = NaN;

    % Track max feasible torque considering only grid points above zero speed
    feasible_torque_values = Te_grid_SI(feasible_mask & Omega_m_grid_SI > 0);
    if isempty(feasible_torque_values)
        max_torque_feasible(k) = 0;
    else
        max_torque_feasible(k) = max(feasible_torque_values);
    end

    % ---- Plot
    fig = figure; hold on; grid on;
    hCF = contourf(Omega_m_grid_SI, Te_grid_SI, Pmech_SI/1e3, 20, 'ShowText', true);
    cb = colorbar;
    ylabel(cb, 'Mechanical power  P_m  (kW)');
    max_kw = max(Pmech_SI(feasible_mask)/1e3, [], 'omitnan');
    if ~isfinite(max_kw) || max_kw <= 0
        max_kw = 1;
    end
    caxis([0, max_kw]);
    xlabel('\omega_m  (mechanical rad/s)');
    ylabel('Electromagnetic torque  T_e  (N\cdotm)');
    title(sprintf('Mechanical Power Contours  P_m  (kW)  |  I_f = %.3f A', If));
    xlim([mechanical_velocity_om_m(1), mechanical_velocity_om_m(end)]);
    ylim([Te_values_Nm(1), Te_values_Nm(end)]);

    % Highlight feasible boundary by overlaying current and voltage limits
    contour(Omega_m_grid_SI, Te_grid_SI, V0_pu_grid, [V0_lim_pu, V0_lim_pu], 'LineColor', 'k', 'LineWidth', 1.2);
    contour(Omega_m_grid_SI, Te_grid_SI, I_mag_pu, [I_lim_pu, I_lim_pu], 'LineColor', [0.3 0.3 0.3], 'LineWidth', 1.2, 'LineStyle', '--');

    % Print summary for this case
    fprintf('If = %.3f A  ->  max feasible torque within grid: %.1f N*m\n', If, max_torque_feasible(k));

    % Save graph if enabled
    if save_graphs && isvalid(fig)
        filename = sprintf('part_2_plot_Pm_contours_If_%0.3fA.png', If);
        saveas(fig, fullfile(outputs_directory, filename));
    end
end

fprintf('\nPart 2 contour plots generated.\n');
