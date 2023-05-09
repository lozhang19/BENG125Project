function dFSHdt = FSH_equation_sensitivity_gamaF(t, FSH, E2, Ih, gammaF)
    % Parameters
    dF = 0.4;
    alphaF = 0.3;
    rF= 1.5;

    % Equation
    dFSHdt = rF * (1 + alphaF * E2(t)) * (1 - (Ih(t)/(gammaF + Ih(t)))) - dF * FSH;
end