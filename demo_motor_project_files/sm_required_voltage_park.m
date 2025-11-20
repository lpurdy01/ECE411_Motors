function [V, Vq, Vd, V0, electrical_angular_velocity, Zs, E] = sm_required_voltage_park(mechanical_velocity_om_m, If, Iq, Id, P)
% sm_required_voltage
% KVL in qd:  V_qd = E_qd + (Rs + j*ωe*Ls)*I_qd
% Lipo–Novotny: q = real{·}, d = −imag{·}

% electrical speed (vector ok)
electrical_angular_velocity = (P.p/2).*mechanical_velocity_om_m;

% stator impedance and internal EMF phasor
Zs = P.Rs + 1i*electrical_angular_velocity.*P.Ls;
E  = (electrical_angular_velocity .* P.Lm .* If) + 1i*0;   % Eq = ωe*Lm*If, Ed = 0

% current phasor (Id is Lipo–Novotny → −imag)
I = Iq + 1i*(-Id);

% required terminal voltage
V = E + Zs.*I;

% split into components (LN peak)
Vq = real(V);
Vd = -imag(V);
V0 = abs(V);
end
