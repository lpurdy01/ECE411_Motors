function P = machine_params()
    % Defines the machine parameters for the synchronous machine.
    % Outputs a struct P with the parameters.
    % --- machine paramters ---
    P.Lm  = 17.5e-3;     % magnetizig induct H
    P.Ls  = 0.6e-3;    % stator induct (q=d) [H]
    P.Rs  = 0.0155;      % stator resist ohm
    P.Rf  = 26.5;      % field resist [ohm]
    P.p   = 8;          % poles
    P.wb  = 2*pi*267;    % elec base w rad/s
    P.wmb = 419;        % mech base w rad/s

    % --- inverter / field ctrl ---
    P.Vdc    = 550;      % dc bus V
    P.I0_max = 300;   % max ac curr (phase peak) A
    P.If_max = 7.0;      % max field curr A dc

    % --- derived ---
    P.Xs_at_elec_base = P.wb*P.Ls;                  % sync react ohm
    P.Xm_at_elec_base = P.wb*P.Lm;                % mag react ohm
    P.f_e_at_elect_base = P.wb/(2*pi);             % elec base freq hz
    P.n_sync_rpm_at_elect_base = P.wmb*60/(2*pi);  % mech base rpm

end