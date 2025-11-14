%% clear
clear; close all; clc;

%% fundamental freq and time vec
f = 60;
we = 2*pi*f;
T = 1/f;
cycles = 3;
npts = 1000;
t = linspace(0, cycles*T, npts);

% amplitude and phase
X_0   = 10;
theta = -90*(pi/180);

% harmonic settings
h           = 1;        % harmonic order (e.g., 3rd)
harm_ratio  = 1/5;      % amplitude ratio relative to fundamental
X_h         = X_0*harm_ratio;

%% abc signals (pure fundamental)
x_a = X_0*cos(we*t + theta);
x_b = X_0*cos(we*t + theta - 2*pi/3);
x_c = X_0*cos(we*t + theta - 4*pi/3);

%% abc signals with harmonic added
% For h = 3, the phase shifts wrap by 2*pi -> the 3rd harmonic is zero-sequence.
x_a_harm = X_0*cos(we*t + theta)                      + X_h*cos(h*we*t + theta);
x_b_harm = X_0*cos(we*t + theta - 2*pi/3)             + X_h*cos(h*we*t + theta - h*2*pi/3);
x_c_harm = X_0*cos(we*t + theta - 4*pi/3)             + X_h*cos(h*we*t + theta - h*4*pi/3);

%% plots: fundamental vs with harmonic
figure
plot(t, x_a, t, x_b, t, x_c, 'LineWidth', 1.1)
title('abc signals (fundamental only)')
xlabel('Time [s]'); ylabel('Amplitude')
legend('x_a','x_b','x_c'); grid on

figure
plot(t, x_a_harm, t, x_b_harm, t, x_c_harm, 'LineWidth', 1.1)
title(sprintf('abc signals with %d^{th} harmonic (%.1f%% of fundamental)', h, 100*harm_ratio))
xlabel('Time [s]'); ylabel('Amplitude')
legend('x_{a,harm}','x_{b,harm}','x_{c,harm}'); grid on

%% Clarke transform (complex space vector form)
j = 1i;
a = exp(j*2*pi/3);

% Use the harmonic-injected signals for the transforms
xa = x_a_harm; xb = x_b_harm; xc = x_c_harm;

x_alphabeta = (2/3)*(xa + a*xb + a^2*xc);   % zero-sequence (e.g., 3rd harmonic) cancels here
x_alpha = real(x_alphabeta);
x_beta  = -imag(x_alphabeta);

% Plot the alpha-beta signals
figure
plot(t, x_alpha, t, x_beta, 'LineWidth', 1.1)
title('\alpha-\beta Signals (Clarke)')
xlabel('Time [s]'); ylabel('Amplitude')
legend('x_{\alpha}', 'x_{\beta}'); grid on

% Space vector trajectory
figure
plot3(x_alpha, -x_beta, t, 'LineWidth', 1.0) % NOTICE -x_beta
title('Clarke Space Vector')
xlabel('\alpha'); ylabel('-\beta'); zlabel('Time [s]')
grid on

%% Park transform (rotate by fundamental angle)
theta_ref = we*t;
x_qd = x_alphabeta .* exp(-j*theta_ref);

x_q = real(x_qd);
x_d = -imag(x_qd);    % negative imaginary (convention)

figure
plot(t, x_q, t, x_d, 'LineWidth', 1.1)
title('q–d Components in \omega_e t Reference Frame (Park)')
xlabel('Time [s]'); ylabel('Amplitude')
legend('x_q','x_d'); grid on

theta_at_end = (180/pi)*atan2(-x_d(1), x_q(1))  % NOTICE:
% 1) -d component
% 2) atan2 function (retains -/- sign)
