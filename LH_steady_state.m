function dLHdt_ss = LH_steady_state(LH, E2_value, P4_value, rL)

    % Parameters 
    dL = 0.4;
    alphaL = 0.3;
    alphaL2 = 0.3;

    % Steady-state equation
    dLHdt_ss = rL * (1 + alphaL * E2_value - alphaL2 * P4_value) - dL * LH;
end
