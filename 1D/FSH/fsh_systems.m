function dydt = fsh_systems(t, y, p)
    % Parameters
    rF = 1.5;
    dF = p(1);
    alphaF = 0.3;
    gammaF = 0.1;

    e0 = 37.4;
    e1 = 150;
    e2 = 12;
    e3 = 115;
    e4 = 18;
    
    i0 = 0.4;
    i1 = 4;
    i2 = 6;
    
    E2 = @(t) e0 + e1 * exp(-(t-14).^2/e2) + e3 * exp(-(t-23).^2/e4) + e1 * exp(-(t-45).^2/e2) + e3 * exp(-(t-54).^2/e4);
    Ih = @(t) i0 + i1 * exp(-(t-22)^2/i2) + i2 * exp(-(t-53)^2/i2);
    
    % FSH equation
    dydt = rF * (1 + alphaF * E2(t)) * (1 - (Ih(t) / (gammaF + Ih(t)))) - dF * y;
end