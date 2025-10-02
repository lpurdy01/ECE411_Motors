% File: hw2_pole_check_from_den.m
% Purpose: Compute poles from the transfer-function denominator and
%          cross-check with the A-matrix eigenvalues.

clear; clc;

%% Motor params (rated flux)
R = 1.43;              % [ohm]
L = 2.38e-3;           % [H]
K = 0.0991;            % [N·m/A] = [V·s/rad]
J = 5.508e-5;          % [kg·m^2]  (converted from 0.0078 oz·in·s^2)
B = 0;                 % [N·m·s/rad]  (set >0 if you want to test damping)

%% Denominator of G(s) = omega(s)/V(s)
% General coupled PMDC (open-loop) has:
%   s^2 + (R/L + B/J) s + (RB/(LJ) + K^2/(LJ))
a2 = 1;
a1 = R/L + B/J;
a0 = (R*B)/(L*J) + (K^2)/(L*J);
den = [a2 a1 a0];

% Poles from the TF denominator (roots)
p_tf = roots(den);

% Cross-check: eigenvalues of the state matrix
A = [-R/L, -K/L;
      K/J, -B/J];
p_eig = eig(A);

%% Report
fprintf('Denominator coefficients: [1  %.6f  %.6f]\n', a1, a0);
fprintf('Poles from TF denominator (roots):   %10.4f  %10.4f  [1/s]\n', p_tf);
fprintf('Poles from state matrix (eig(A)):     %10.4f  %10.4f  [1/s]\n', p_eig);

% Optional: show time constants and damping ratio (for intuition)
tau_a = L/R;                     % electrical time constant
tau_m = (J*R)/(K^2);             % electromechanical time constant (B=0 case)
zeta  = (a1)/(2*sqrt(a0));       % damping ratio of the 2nd-order denominator
wn    = sqrt(a0);                % natural frequency

fprintf('\nTau_a = %.4e s,  Tau_m = %.4e s,  zeta = %.3f,  wn = %.1f rad/s\n', ...
        tau_a, tau_m, zeta, wn);
%% ================== Part 2(c): K = 300% of nominal ==================
K3 = 3*K;                     % 300% flux => K scaled by 3

% Denominator of G(s) = ω(s)/V(s) with possible viscous B
a2_3 = 1;
a1_3 = R/L + B/J;
a0_3 = (R*B)/(L*J) + (K3^2)/(L*J);
den3 = [a2_3 a1_3 a0_3];

% Poles from TF denominator and from A-matrix
p_tf_3  = roots(den3);
A3      = [-R/L, -K3/L;
            K3/J, -B/J];
p_eig_3 = eig(A3);

% Damping ratio and natural frequency (2nd-order denominator)
wn3   = sqrt(a0_3);
zeta3 = a1_3/(2*wn3);

fprintf('\n=== Part 2(c): Poles with K = 300%% of nominal (K = %.4f) ===\n', K3);
fprintf('Denominator: s^2 + (%.6f)s + (%.6f)\n', a1_3, a0_3);
fprintf('Poles from TF denominator (roots):   %10.4f %10.4f [1/s]\n', p_tf_3);
fprintf('Poles from state matrix (eig(A)):     %10.4f %10.4f [1/s]\n', p_eig_3);
fprintf('wn = %.2f rad/s,  zeta = %.3f  => %s\n', wn3, zeta3, ...
    ternary(zeta3<1,'UNDERDAMPED (oscillatory)','NON-oscillatory'));

% --- tiny local helper for printing ---
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

%% ================== Part 2(d,e): pole plot + damping metrics ==================
% Two cases: nominal K and 3x nominal
K1  = K;        % 100%
K3  = 3*K;      % 300%

% Denominators: s^2 + (R/L + B/J)s + (RB/(LJ) + K^2/(LJ))
den1 = [1,  R/L + B/J,  (R*B)/(L*J) + (K1^2)/(L*J)];
den3 = [1,  R/L + B/J,  (R*B)/(L*J) + (K3^2)/(L*J)];

p1 = roots(den1);
p3 = roots(den3);

%% ---- Damping & classification (robust) ----
wn1   = sqrt(den1(3));  zeta1 = den1(2)/(2*wn1);
wn3   = sqrt(den3(3));  zeta3 = den3(2)/(2*wn3);

fprintf('\n=== 2(e) Damping & classification ===\n');
fprintf('K=100%%:  zeta = %.3f,  wn = %.1f rad/s  -> %s\n', zeta1, wn1, classify_damping(zeta1));
fprintf('K=300%%:  zeta = %.3f,  wn = %.1f rad/s  -> %s\n', zeta3, wn3, classify_damping(zeta3));

% ---- helper at end of file (after all code) ----
function s = classify_damping(z)
    tol = 1e-3;  % tolerance for "critical"
    if z > 1 + tol
        s = 'Overdamped';
    elseif z < 1 - tol
        s = 'Underdamped';
    else
        s = 'Critically damped';
    end
end

% ---------- Plot Re–Im pole map ----------
figure('Color','w','Name','Poles at K=100% and K=300%'); hold on; grid on;
% Axes
xline(0,'k:'); yline(0,'k:');

% Plot poles
plot(real(p1), imag(p1), 'o', 'MarkerSize',8, 'MarkerFaceColor',[0.10 0.70 0.10], 'MarkerEdgeColor','k'); % green
plot(real(p3), imag(p3), 's', 'MarkerSize',8, 'MarkerFaceColor',[0.10 0.30 0.90], 'MarkerEdgeColor','k'); % blue

% Labels
for k = 1:numel(p1)
    text(real(p1(k))+10, imag(p1(k))+10, 'K=100%', 'Color',[0.10 0.70 0.10]);
end
for k = 1:numel(p3)
    text(real(p3(k))+10, imag(p3(k))+10, 'K=300%', 'Color',[0.10 0.30 0.90]);
end

xlabel('Real\{p\}  [1/s]');
ylabel('Imag\{p\}  [1/s]');
title('Eigenvalues (poles) on the complex plane: K=100% (green) vs K=300% (blue)');
axis equal;  % circles look like circles
% Set a reasonable window around both cases
xr = [min([real(p1);real(p3)])-100, max([real(p1);real(p3)])+100];
yr = [min([imag(p1);imag(p3)])-100, max([imag(p1);imag(p3)])+100];
xlim(xr); ylim(yr);

legend({'K=100% poles','K=300% poles'},'Location','best');
saveas(gcf,'hw2_poles_xy.png');


%% ---------- Root-locus style pole map (centered, square) ----------
figure('Color','w','Name','Poles: K=100% (green) vs K=300% (blue)');
hold on; grid on; box on;

% Collect poles
P = [p1(:); p3(:)];
xr = max(abs(real(P))) * 1.15;
yr = max(abs(imag(P))) * 1.15;
L  = max(xr, yr);                      % symmetric limits

% Axes centered at origin
plot([-L L],[0 0],'k-','LineWidth',1.0,'HandleVisibility','off'); % real axis
plot([0 0],[-L L],'k-','LineWidth',1.0,'HandleVisibility','off'); % imag axis
set(gca,'XAxisLocation','origin','YAxisLocation','origin');

% Plot poles with nice markers
plot(real(p1), imag(p1), 'o', 'MarkerSize',9, 'MarkerFaceColor',[0.10 0.70 0.10], 'MarkerEdgeColor','k', ...
     'DisplayName','K = 100% poles');
plot(real(p3), imag(p3), 's', 'MarkerSize',9, 'MarkerFaceColor',[0.10 0.30 0.90], 'MarkerEdgeColor','k', ...
     'DisplayName','K = 300% poles');

% Annotations (optional)
text(real(p1)+15, imag(p1)+15, "100%", 'Color',[0.10 0.70 0.10]);
text(real(p3)+15, imag(p3)+15, "300%", 'Color',[0.10 0.30 0.90]);

% Tidy axes
axis equal; axis([-L L -L L]);      % square & symmetric
xlabel('Real  [1/s]');
ylabel('Imag  [1/s]');
title('d: poles on the complex plane root locust style');
legend('Location','northeast');
