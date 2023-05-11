function dFSHdt_ss = FSH_steady_state(FSH, E2, Ih, rF, dF)
    % Parameters
    alphaF = 0.3;
    gammaF = 0.1;

    % Steady-state equation
    dFSHdt_ss = rF * (1 + alphaF * E2) * (1 - (Ih / (gammaF + Ih))) - dF * FSH;
end