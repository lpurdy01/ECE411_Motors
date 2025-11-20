function B = pu_bases(VB_LN_rms, PB_3ph, wB_elec, p)
%PU_BASES Returns a struct of base quantities for per-unit analysis.
%   VB_LN_rms : rated line-to-neutral RMS voltage [V]
%   PB_3ph    : rated three-phase power [W]
%   wB_elec   : electrical base angular frequency [rad/s]
%   p         : number of poles (used for reporting mechanical base speed)
%
% Primary bases (LN peak voltage consistent with Park frame usage)
B.VB  = sqrt(2) * VB_LN_rms;   % line-neutral peak base voltage [V]
B.PB  = PB_3ph;                % three-phase base power [W]
B.wB  = wB_elec;               % electrical base angular speed [rad/s]

% Secondary bases derived from the above
B.IB  = (2 * B.PB) / (3 * B.VB);  % phase current peak base [A]
B.ZB  = B.VB / B.IB;              % impedance base [Ohm]
B.LB  = B.ZB / B.wB;              % inductance base [H]
B.wmB = 2 * B.wB / p;             % mechanical base angular speed [rad/s]
B.TB  = (p / 2) * (B.PB / B.wB);  % torque base [N*m]
end
