function dFSHdt = FSH_system(t, x, rF, E2_value, Ih_value)
    FSH = x(1);

    % Parameters
    dF = 0.4;
    alphaF = 0.3;
    gammaF = 0.1;

    % Equation
    dFSHdt(1,1) = rF * (1 + alphaF * E2_value) * (1 - (Ih_value/(gammaF + Ih_value))) - dF * FSH;
end