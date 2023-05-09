function dLHdt = LH_equation_sensitivity_alphaL(t, LH, E2, P4, dL)
    % Parameters
    alphaL = 0.3;
    rL = 1.5;
    alphaL2 = 0.3;

    % Equation
    dLHdt = rL * (1 + alphaL * E2(t)-alphaL2* P4(t)) - dL * LH;
end