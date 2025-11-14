%% prob 5 c
% find base qnts: PB, VB, IB, ZB, TB, wB, wmB
% then perunit: rs_pu, Xs_pu, Ls_pu
% assumptons: LN peak bases, 3ph power base, 60hz, p=4, lipo-novotny normlization

clear; clc;

%% givens (rated data)
V_LN_RMS = 277;       % line-neutral rms volt [V]
P_rated  = 20e3;      % 3phase rated power [W]
f        = 60;    % elec freq hz
p        = 4;         % poles
rs       = 0.180;     % stator resist [ohm]
Ls       = 26e-3;     % sync induct [H]

%% base values (primary)
VB = sqrt(2)*V_LN_RMS;    % volt base (LN peak)
PB = P_rated;  % power base (3ph)
wB = 2*pi*f;              % electrical base speed rad/s

%% base values (secondary)
%  IB = (2/3)*PB/VB
IB = (2/3)*PB/VB;         % current base (LN peak) A
ZB = VB/IB;               % imped base ohm
LB = ZB/wB;               % induct base H
wmB = (2/p)*wB;           % mech base speed rad/s
TB  = PB/wmB;             % torque base N*m

%% machine params (per unit @ 60hz)
Xs    = wB*Ls;            % sync react ohm
rs_pu = rs/ZB;            % pu stator res
Xs_pu = Xs/ZB;    % pu sync react
Ls_pu = Ls/LB;            % pu induct (shd ~= Xs_pu)

%% print stuyle results
fprintf('p 5c perunit params\n');
fprintf('Bases:\n');
fprintf('  PB  = %.3f kW\n', PB/1e3);
fprintf('  VB  = %.3f Vpeak (LN) \n', VB);
fprintf('  IB  = %.3f Apeak (LN)\n', IB);
fprintf('  ZB   = %.6f ohms\n', ZB);
fprintf('  LB  = %.6e H\n', LB);
fprintf('  TB   = %.3f Nm\n', TB);
fprintf('  wB   = %.3f rad/s elec\n', wB);
fprintf('  wmB  = %.3f rad/s mech\n', wmB);

fprintf('\per-unit prams 60hz:\n');
fprintf('  rs_pu = %.6f pu\n', rs_pu);
fprintf('  Xs_pu = %.6f pu \n', Xs_pu);
fprintf('  Ls_pu = %.6f pu  (chk ~ Xs_pu)\n', Ls_pu);
