% New function for LH_derivative
function dLHdt = LH_derivative(LH, t, E2, P4)
    % Parameters
    rL = 1.5;
    dL = 0.4;
    alphaL = 0.3;
    alphaL2 = 0.3;

    % Equation
    dLHdt = rL * (1 + alphaL * E2(t) - alphaL2 * P4(t)) - dL * LH;
end