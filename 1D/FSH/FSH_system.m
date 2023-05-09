function dFSHdt = FSH_system(x, p)
    FSH = x;
    rF = p(1);
    E2_value = p(2);
    Ih_value = p(3);

    % Parameters
    dF = 0.4;
    alphaF = 0.3;
    gammaF = 0.1;

    % Equation
    dFSHdt = rF * (1 + alphaF * E2_value) * (1 - (Ih_value/(gammaF + Ih_value))) - dF * FSH;
end