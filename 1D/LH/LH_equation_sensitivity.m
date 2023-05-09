function dLHdt = LH_equation_sensitivity(t, LH, E2, P4, rL)
    % Parameters
    dL = 0.4;
    alphaL = 0.3;
    alphaL2 = 0.3;

    % Equation
    dLHdt = rL * (1 + alphaL * E2(t)-alphaL2* P4(t)) - dL * LH;
end