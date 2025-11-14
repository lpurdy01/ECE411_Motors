%% Clear cloose
clear; clc; close all;

%% machine + sim params
Lambda0 = 10;              % flux constant [V*s/rad]
we = -1.0;                 % electrical speed rad/s (neg = clockwise)

% time stuff
T = 2*pi/abs(we);          % period [s]
npts = 1000;               % points per period
t = linspace(0, T, npts);  % time vec

%% rotor flux in alpha-beta
j = 1i;                                        % imaginary unit
lambda_alpha_beta = Lambda0 * exp(j * we * t); % space vector

% split into components
lambda_alpha = real(lambda_alpha_beta);        % alpha comp
lambda_beta  = -imag(lambda_alpha_beta);       % beta comp (neg imag)

%% stator back emf in _ab 
e_alpha_beta = j * we * lambda_alpha_beta;

% split components
e_alpha = real(e_alpha_beta);
e_beta  = -imag(e_alpha_beta);

%% plots — rotor flux + back emf in ab
figure
plot3(lambda_alpha, lambda_beta, t); hold on
plot3(e_alpha, e_beta, t, '--');
grid on
xlabel('alpha (real)')
ylabel('beta (neg imag)')
zlabel('time (s)')
title('rotor flux and back emf in alpha-beta')
legend('lambda_ab','e_ab')
view(30,30)

%% next section

%% 
% part (b) — qd (park) rotor flux and back emf
% lipo–novotny: q || alpha at t=0; q = real{·}, d = -imag{·}
% park ref angle: theta_ref = -we*t
% qd should be constant in time

% park transform complex form): f_qd = f_ab * exp(-j*theta_ref)
j = 1i;
theta_ref = we .* t;

lambda_qd = lambda_alpha_beta .* exp(-j .* theta_ref);
e_qd      = e_alpha_beta      .* exp(-j .* theta_ref);

% pull out q,d (using class mapping)
lambda_q = real(lambda_qd);
lambda_d = -imag(lambda_qd);

e_q = real(e_qd);
e_d = -imag(e_qd);

% xy plot with q on +x and -d on +y
ptL = [mean(lambda_q), -mean(lambda_d)];
ptE = [mean(e_q),      -mean(e_d)];

figure; hold on; grid on; box on;
title('part (b): qd rotor flux and back emf (theta = we*t)')
xlabel('q comp'); ylabel('-d comp'); axis equal

mx = max(abs([ptL(1), ptE(1), 0]));
my = max(abs([ptL(2), ptE(2), 0]));
sx = max(1, 1.3*mx);
sy = max(1, 1.3*my);
xlim([-sx sx]); ylim([-sy sy]);

plot([-sx sx],[0 0],'k:')
plot([0 0],[-sy sy],'k:')

scatter(ptL(1), ptL(2), 60, 'filled')
scatter(ptE(1), ptE(2), 60, '^', 'filled')

legend('Lambda_qd','E_qd','Location','southeast')

% quick console print
fprintf('part (b): Lambda q=%.6g, d=%.6g (plot -d=%.6g)\n', mean(lambda_q), mean(lambda_d), -mean(lambda_d));
fprintf('part (b): E      q=%.6g, d=%.6g (plot -d=%.6g)\n', mean(e_q), mean(e_d), -mean(e_d));
