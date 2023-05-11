function dFSHdt = FSH_equation_sensitivity_dF(t, FSH, E2, Ih, dF)
    % Parameters
    rF = 0.4;
    alphaF = 0.3;
    gammaF = 0.1;

    % Equation
    dFSHdt = rF * (1 + alphaF * E2(t)) * (1 - (Ih(t)/(gammaF + Ih(t)))) - dF * FSH;
end
