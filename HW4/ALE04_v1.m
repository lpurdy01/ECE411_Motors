%% clear
clear; close; clc;

%% fundamentla freq and time vec

f = 60;
we = 2*pi*f;
T = 1/f;
cycles = 3;
npts = 1000;
t = linspace(0,cycles*T,npts);

% amp and phase
X_0 = 10;
theta = 120*(pi/180);

%% abc signals

x_a = X_0*cos(we*t+theta) + 0;
x_b = X_0*cos(we*t+theta-2*pi/3) + 0;
x_c = X_0*cos(we*t+theta-4*pi/3) + 0;

figure
plot(t,x_a,t,x_b,t,x_c)
title('abc signals')
xlabel('Time [s]')
ylabel('Amplitude')
legend('x_a','x_b','x_c')

%% clark transform

j = 1i;
a = exp(j*2*pi/3);

x_alphabeta = (2/3)*(x_a + a*x_b + x_c*a^2);
x_alpha = real(x_alphabeta);
x_beta = -imag(x_alphabeta);

% Plot the alpha-beta signals
figure
plot(t, x_alpha, t, x_beta)
title('Alpha-Beta Signals')
xlabel('Time [s]')
ylabel('Amplitude')
legend('x_{alpha}', 'x_{beta}')

figure
plot3(x_alpha,-x_beta,t) % NOTICE -x_beta
title('Clarke space vector')
xlabel('Alpha')
ylabel('-Beta')
zlabel('Time [s]')
grid on

%% Park transform

theta_ref = we*t;
x_qd = x_alphabeta.*exp(-j*theta_ref);

x_q = real(x_qd);
x_d = -imag(x_qd); % negative imaginary

figure
plot(t,x_q,t,x_d)
title('qd Components in we*t ref frame')
xlabel('Time [s]')
ylabel('Amplitude')
legend('x_q','x_d')

theta_at_end = (180/pi)*atan2(-x_d(1),x_q(1)) % NOTICE:
% 1) - d component
% 2) atan2 function (retains -/- sign)