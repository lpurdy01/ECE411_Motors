%% Clear all close all Problem 4
% notes
% Convention: Lipo–Novotny (q-axis ∥ α at θ=0)
% Reverse Park (T2^{-1}): [ i_α ; i_β ] = [ cosθ  sinθ ; -sinθ  cosθ ] [ i_q ; i_d ]
% Reverse Clarke (T1^{-1}, balanced so f0=0):
%   [ i_a ]   [ 1      0     ] [ i_α ]
%   [ i_b ] = [ -1/2  -√3/2 ] [ i_β ]
%   [ i_c ]   [ -1/2   √3/2 ] [     ]
% (If you keep zero-sequence explicitly: multiply [i_α i_β i_0]^T by
%   [ 1  0  1 ; -1/2 -√3/2 1 ; -1/2 √3/2 1 ], with i_0 = 0 for balanced 3Φ.)

clear; close all; clc;

%% Common timing
f_e = 60;                    % electrical frequency [Hz]
T   = 1/f_e;                 % given
t_end = 3*T;                 % 3 electrical periods
npts  = 2000;
t = linspace(0, t_end, npts);

%% -----------------------------
%  (a) q- and d-axis currents vs time (as given)
tau = 0.005;                 % seconds
i_q = 1.0 * (1 - exp(-t./tau));
i_d = 2.0 * (1 - exp(-t./tau));

figure; hold on; grid on; box on;
plot(t, i_q, 'LineWidth', 2);
plot(t, i_d, 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Current [A]');
title('4(a): q- and d-axis Currents vs Time (Park frame)');
legend('i_q(t) = 1\,(1 - e^{-t/\tau})', 'i_d(t) = 2\,(1 - e^{-t/\tau})', 'Location','southeast');

%% -----------------------------
%  (b) α,β currents vs time (reverse Park transform, θ = 2π·60·t)
theta = 2*pi*f_e*t;          % Park angle (per prompt)
c = cos(theta); s = sin(theta);

% Reverse Park: T2^{-1}(θ)
i_alpha =  c .* i_q + s .* i_d;
i_beta  = -s .* i_q + c .* i_d;

figure; hold on; grid on; box on;
plot(t, i_alpha, 'LineWidth', 2);
plot(t, i_beta,  'LineWidth', 2);
xlabel('Time [s]'); ylabel('Current [A]');
title('4(b): \alpha,\beta Currents vs Time (Reverse Park, \theta = 2\pi 60 t)');
legend('i_\alpha(t)', 'i_\beta(t)', 'Location','southeast');

%% -----------------------------
%  (c) a,b,c currents vs time (reverse Clarke of part b)
% Balanced, so i0 = 0
i_a = 1      * i_alpha + 0           * i_beta;
i_b = (-1/2) * i_alpha + (-sqrt(3)/2)* i_beta;
i_c = (-1/2) * i_alpha + ( sqrt(3)/2)* i_beta;

figure; hold on; grid on; box on;
plot(t, i_a, 'LineWidth', 2);
plot(t, i_b, 'LineWidth', 2);
plot(t, i_c, 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Current [A]');
title('4(c): a,b,c Currents vs Time (Reverse Clarke of \alpha\beta)');
legend('i_a(t)','i_b(t)','i_c(t)', 'Location','southeast');

%% Sanity notes:
% • (a) shows first-order rises to 1 A (q) and 2 A (d).
% • (b) αβ are those same transients rotated by θ, so they appear sinusoidal
%   envelopes modulated by the exponential rise.
% • (c) abc are three sinusoids 120° apart with the same exponential envelopes.
