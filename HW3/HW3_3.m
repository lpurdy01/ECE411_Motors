%% Clear cloose
clear; clc; close all;

%% Machine and simulation parameters
% Fundamental electrical speed and flux constant
Lambda0 = 10;             % V·s/rad (flux constant)
we = -1.0;                % electrical speed [rad/s] (negative = clockwise)

% Derived time-domain settings
T = 2*pi / abs(we);       % One electrical period [s]
npts = 1000;             % Number of time points per period
t = linspace(0, T, npts); % Time vector

%% Rotor flux in alpha-beta frame
% Using complex space vector notation with Lipo-Novotny convention
j = 1i;                   % imaginary unit
lambda_alpha_beta = Lambda0 * exp(j * we * t);   % Rotor flux space vector

% Separate real and imaginary components (alpha and beta axes)
lambda_alpha = real(lambda_alpha_beta);         % alpha-axis component
lambda_beta  = -imag(lambda_alpha_beta);        % beta-axis component (negative sign per class convention)

%% Stator back-EMF in alpha-beta frame
% e = dλ/dt = j*we * λab
e_alpha_beta = j * we * lambda_alpha_beta;

% Separate components
e_alpha = real(e_alpha_beta);
e_beta  = -imag(e_alpha_beta);                  % Again, negative sign convention

%% Plot results - Rotor Flux and Back-EMF
figure
plot3(lambda_alpha, lambda_beta, t, 'LineWidth', 2); hold on
plot3(e_alpha, e_beta, t, '--', 'LineWidth', 2);

grid on
xlabel('\alpha-axis (real)'); 
ylabel('-\beta-axis (imag)'); 
zlabel('Time [s]');
title('Rotor Flux and Back-EMF Space Vectors in \alpha\beta Frame');
legend('\lambda_{\alpha\beta}', 'e_{\alpha\beta}');
view(30,30)

%% Notes:
% - Rotor flux rotates clockwise (negative we), so vector rotates opposite the usual direction.
% - q-axis aligns with alpha at t = 0 per Lipo-Novotny definition.
% - Beta components use negative imaginary parts to match course plotting convention.

%% next section

%% ============================================================
%  Part (b) — q-d (Park) rotor flux Λ_qd and back-EMF E_qd
%  Convention: Lipo–Novotny (q ∥ α at θ=0; q = Real{·}, d = −Imag{·})
%  Park reference angle per prompt: theta_ref = -we * t
%  q-d components should be constant in time in this frame.
% =============================================================

% Safety checks (optional—comment out if variables are guaranteed)
assert(exist('lambda_alpha_beta','var')==1, 'lambda_alpha_beta not found');
assert(exist('e_alpha_beta','var')==1,     'e_alpha_beta not found');
assert(exist('we','var')==1 && exist('t','var')==1, 'we or t missing');

% Park transform in complex form: f_qd = f_ab * exp(-j*theta_ref)
j = 1i;
theta_ref = we .* t;  % <- as required by the problem statement

lambda_qd = lambda_alpha_beta .* exp(-j .* theta_ref);
e_qd      = e_alpha_beta      .* exp(-j .* theta_ref);

% Extract q and d components (Lipo–Novotny mapping)
lambda_q = real(lambda_qd);
lambda_d = -imag(lambda_qd);

e_q = real(e_qd);
e_d = -imag(e_qd);

% ---- XY plot (q on +X, -d on +Y), both vectors as arrows/points ----
% Compute constant points
ptL = [mean(lambda_q), -mean(lambda_d)];   % [ +q ,  -d ]
ptE = [mean(e_q),      -mean(e_d)];

% Figure with nice formatting
figure; hold on; grid on; box on;
title('Part (b): q–d Rotor Flux and Back-EMF (Park frame, \theta = \omega_e t)');
xlabel('+q component'); ylabel('−d component'); axis equal;

% Set sensible axes around the points
mx = max(abs([ptL(1), ptE(1), 0]));
my = max(abs([ptL(2), ptE(2), 0]));
sx = max(1, 1.3*mx);   % add ~30% margin, min span ~1
sy = max(1, 1.3*my);
xlim([-sx sx]); ylim([-sy sy]);

% Draw zero axes
plot([-sx sx],[0 0],'k:','LineWidth',1);
plot([0 0],[-sy sy],'k:','LineWidth',1);


s1 = scatter(ptL(1),ptL(2),70,[0 0.45 0.74],'filled');
s2 = scatter(ptE(1),ptE(2),70,[0.85 0.33 0.10],'filled','^');


legend([s1 s2], '\Lambda_{qd}', 'E_{qd}', 'Location','southeast');

% Console summary (nice to keep)
fprintf('Part (b) constants:\n  Lambda: q=%.6g, d=%.6g  (plot uses -d=%.6g)\n  EMF:    q=%.6g, d=%.6g  (plot uses -d=%.6g)\n', ...
        mean(lambda_q), mean(lambda_d), -mean(lambda_d), mean(e_q), mean(e_d), -mean(e_d));
