function [V_pu, Vq_pu, Vd_pu, V0_pu, w_e, Zs_pu, E_pu] = ...
    sm_required_voltage_park_pu(om_m, If, Iq, Id, P, B)
%SM_REQUIRED_VOLTAGE_PARK_PU Per-unit Park-frame voltage calculation.
%   om_m : mechanical angular speed [rad/s] (vector ok)
%   If   : field current [A dc]
%   Iq   : q-axis stator current [A peak]
%   Id   : d-axis stator current [A peak]
%   P    : machine parameter struct from machine_params()
%   B    : base quantities struct from pu_bases()
%
% Returns:
%   V_pu  : complex phase voltage in per-unit
%   Vq_pu : q-axis component (LN peak) in per-unit
%   Vd_pu : d-axis component (LN peak) in per-unit
%   V0_pu : magnitude in per-unit
%   w_e   : electrical angular speed [rad/s]
%   Zs_pu : stator impedance in per-unit
%   E_pu  : internal EMF in per-unit

% electrical speed (vector ok)
w_e  = (P.p/2) .* om_m;
w_pu = w_e ./ B.wB;

% base impedance conversions
Rs_pu = P.Rs / B.ZB;
Xs_base_pu = (B.wB * P.Ls) / B.ZB;

% speed-scaled reactances (per-unit)
Zs_pu = Rs_pu + 1i * (w_pu .* Xs_base_pu);
E_pu  = (w_e .* P.Lm .* If) / B.VB + 1i*0;

% stator current in per-unit (Lipo–Novotny: d = -imag)
I_pu = (Iq + 1i * (-Id)) / B.IB;

% required stator voltage in per-unit
V_pu = E_pu + Zs_pu .* I_pu;

% split into components and magnitude
Vq_pu = real(V_pu);
Vd_pu = -imag(V_pu);
V0_pu = abs(V_pu);
end
