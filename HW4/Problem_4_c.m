%% Clear Close section
clc; clear; close;

% Setup variables
Ld=0.9;E0=1.25;V0=1;we=1;
Lq=0.01:0.001:0.9;
Te_field=E0*V0./(we*Ld);
Te_rel=V0^2*(Ld-Lq)./(2*we*Ld.*Lq);
min_total_field = min(abs(Te_rel-Te_field));
idx=find(abs(Te_rel-Te_field)==min_total_field);
Lq_eq=Lq(idx);
Sal= Ld/Lq_eq;

% Output the results
fprintf('The value of Lq that makes the maximum reluctance torque equal to the maximum field weakening torque is: %.4f pu\n', Lq_eq);
fprintf('The corresponding saliency ratio (Ld/Lq) is: %.4f\n', Sal);    
