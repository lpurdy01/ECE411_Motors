%% Clear Close section
clc; clear; close;

% Setup variables
Ld=0.9;Lq=0.6;E0=1.25;V0=1;we=1;
% Setup range to calc over
delta=linspace(-pi,pi,1000);
% Calculate torque
Te_field=E0*V0./(we*Ld).*sin(delta);
Te_rel=V0^2*(Ld-Lq)./(2*we*Ld*Lq).*sin(2*delta);
Te=Te_field+Te_rel;
Te_field_max=max(Te_field);
Te_rel_max=max(Te_rel);
% plot
figure;
plot(delta,Te_field,delta,Te_rel,delta,Te);

% Label the axes
xlabel('Rotor Angle, \delta (rad)');
ylabel('Electromagnetic Torque, T_e pu');
legend('Field Weakening Torque','Reluctance Torque','Total Torque','Location','Best');
title('Electromagnetic Torque vs Rotor Angle');
grid on;

% Output the max torques
fprintf('Maximum Field Weakening Torque: %.4f pu\n', Te_field_max);
fprintf('Maximum Reluctance Torque: %.4f pu\n', Te_rel_max);
fprintf('Maximum Total Torque: %.4f pu\n', max(Te));
