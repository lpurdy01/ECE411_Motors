%% Clear all close all Problem 4
% notes: lipo-novotny (q || alpha at t=0), reverse park then reverse clarke
clear; close all; clc;

%% common timing
f_e = 60;            % electrical frequency [Hz]
T   = 1/f_e;         % period s
t_end = 3*T;         % 3 periods
npts  = 2000;
t = linspace(0, t_end, npts);

%% ----------------------
% (a) q and d currents vs time (as given)
tau = 0.005;         % seconds
i_q = 1.0 * (1 - exp(-t./tau));
i_d = 2.0 * (1 - exp(-t./tau));

figure; hold on; grid on;
plot(t, i_q, 'LineWidth', 2);
plot(t, i_d, 'LineWidth', 2);
xlabel('time (s)'); ylabel('current (A)');
title('4a: q and d currents vs time (park frame)');
legend('iq(t)','id(t)','Location','southeast');

%% --------------
% (b) alpha,beta currents (reverse park, theta = 2*pi*60*t)
theta = 2*pi*f_e*t;   % park angle
c = cos(theta); s = sin(theta);

% reverse park
i_alpha =  c .* i_q + s .* i_d;
i_beta  = -s .* i_q + c .* i_d;

figure; hold on; grid on;
plot(t, i_alpha, 'LineWidth', 2);
plot(t, i_beta,  'LineWidth', 2);
xlabel('time (s)'); ylabel('current (A)');
title('4b: alpha,beta currents vs time (reverse park)');
legend('i_alpha(t)','i_beta(t)','Location','southeast');

%% -----------------------------
% (c) a,b,c currents (reverse clarke, balanced so i0=0)
i_a = (1)      * i_alpha + (0)             * i_beta;
i_b = (-1/2)   * i_alpha + (-sqrt(3)/2)    * i_beta;
i_c = (-1/2)   * i_alpha + ( sqrt(3)/2)    * i_beta;

figure; hold on; grid on;
plot(t, i_a, 'LineWidth', 2);
plot(t, i_b, 'LineWidth', 2);
plot(t, i_c, 'LineWidth', 2);
xlabel('time (s)'); ylabel('current (A)');
title('4c: a,b,c currents vs time (reverse clarke)');
legend('i_a(t)','i_b(t)','i_c(t)','Location','southeast');

%% quick notes
% (a) first-order rises to 1 A q and 2 A (d)
% b alpha beta are those transients rotated by theta
% (c) three phase sinuoids 120 deg apart with same envelopes
