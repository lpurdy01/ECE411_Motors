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

% --- choose per-unit bases (match earlier assumptions) ---
VB_LN_peak = P.Vdc/sqrt(3);               % inverter linear-modulation ceiling (LN peak)
VB_LN_rms  = VB_LN_peak/sqrt(2);
PB_3ph     = 1.5 * VB_LN_peak * P.I0_max;  % ensures I_base equals I0_max (peak)
B          = pu_bases(VB_LN_rms, PB_3ph, P.wb, P.p);

% handy limit
M0_limit  = 2/sqrt(3);                     % linear region ceiling for M0
k_M0      = 2/sqrt(3);                     % M0 = k_M0 * V0_pu with this base choice

%% --- Plot 3 ------
%
% ██████  ██       ██████  ████████     ██████  
% ██   ██ ██      ██    ██    ██             ██ 
% ██████  ██      ██    ██    ██         █████  
% ██      ██      ██    ██    ██             ██ 
% ██      ███████  ██████     ██        ██████
%                                            

% Speed sweep (mechanical rad/s); 0 → 2*wmb per spec
mechanical_velocity_om_m = linspace(0, 3*P.wmb, 200);  % mechanical rad/s

% Torque sweep (N*m); match Part 2 range
Te_values_Nm = linspace(0, 220, 200);

% Meshgrids and per-unit torque grid
[Omega_m_grid_SI, Te_grid_SI] = meshgrid(mechanical_velocity_om_m, Te_values_Nm);
Te_grid_pu = Te_grid_SI / B.TB;

for k = 1:numel(field_current_cases)
    If = field_current_cases(k);

    % Electrical speeds (SI & pu)
    Omega_e_grid_SI = (P.p/2) .* Omega_m_grid_SI;
    Omega_e_grid_pu = Omega_e_grid_SI / B.wB;

    % Internal EMF in pu (Eq only; Ed=0 here)
    Eq_grid_pu = (Omega_e_grid_SI .* P.Lm .* If) / B.VB;

    % avoid divide-by-zero/undefined
    small_emf_mask = abs(Eq_grid_pu) < 1e-9;
    Eq_grid_pu(small_emf_mask) = NaN;

    % Required Iq from torque definition (Id=0)
    Iq_grid_pu = (Te_grid_pu .* Omega_e_grid_pu) ./ Eq_grid_pu;
    Iq_grid_pu(~isfinite(Iq_grid_pu)) = NaN;

    % Convert currents (pu -> SI) for helper function
    Iq_grid_SI = Iq_grid_pu * B.IB;
    Id_grid_SI = zeros(size(Iq_grid_SI));

    % Flatten, evaluate, reshape
    Iq_vec = Iq_grid_SI(:);
    Id_vec = Id_grid_SI(:);
    Om_vec = Omega_m_grid_SI(:);

    [~, ~, ~, V0_pu_vec] = sm_required_voltage_park_pu( ...
        Om_vec, If, Iq_vec, Id_vec, P, B);

    V0_pu_grid = reshape(V0_pu_vec, size(Omega_m_grid_SI));

    % Modulation index (dimensionless)
    M0_grid = k_M0 * V0_pu_grid;

    % ---- Plot
    fig = figure; hold on; grid on;
    contourf(Omega_m_grid_SI, Te_grid_SI, M0_grid, 20, 'ShowText', true);
    cb = colorbar; ylabel(cb, 'Modulation Index  M_0');

    % same axes ranges as before for comparisons
    xlabel('\omega_m  (mechanical rad/s)');
    ylabel('Electromagnetic torque  T_e  (N\cdotm)');
    title(sprintf('Inverter Modulation Amplitude  M_0  |  I_f = %.3f A', If));
    xlim([mechanical_velocity_om_m(1), mechanical_velocity_om_m(end)]);
    ylim([Te_values_Nm(1), Te_values_Nm(end)]);

    % overlay linear-modulation limit curve
    contour(Omega_m_grid_SI, Te_grid_SI, M0_grid, [M0_limit M0_limit], ...
        'LineColor', 'k', 'LineWidth', 1.2);                     % M0 = 2/sqrt(3)

    % Save per-case plot
    if save_graphs && isvalid(fig)
        filename = sprintf('part_3_plot_M0_contours_If_%0.3fA.png', If);
        saveas(fig, fullfile(outputs_directory, filename));
    end
end

fprintf('\nPart 3 modulation contour plots generated.\n');
