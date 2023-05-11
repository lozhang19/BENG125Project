function dydt = lh_systems(t, y, e1)
    % Parameters
    rL = 1.5;
    dL = 0.4;
    alphaL = 0.3;
    alphaL2 = 0.3;
    
    e0 = 37.4;
    e2 = 12;
    e3 = 115;
    e4 = 18;
    
    p0 = 30;
    p1 = 15;
    p2 = 15;
    
    E2 = @(t) e0 + e1 * exp(-(t-14).^2/e2) + e3 * exp(-(t-23).^2/e4) + e1 * exp(-(t-45).^2/e2) + e3 * exp(-(t-54).^2/e4);
    P4 = @(t) p0 + p1 * exp(-(t-22)^2/p2) + p1 * exp(-(t-53)^2/p2);
    
    % LH equation
    dydt = rL * (1 + alphaL * E2(t) - alphaL2 * P4(t)) - dL * y;
end