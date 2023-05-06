function dydt = hormone_system(t, y, rF, rL, rE, rP, dF, dL, dE, dP, alphaF, alphaL, KmE, KmP, tauH, tauG)
    H = 1 + 0.5 * sin(2 * pi * t / tauH);
    G = 1 + 0.5 * sin(2 * pi * t / tauG);
    
    dydt = [rF * H * (1 + alphaF * y(3)) - dF * y(1);
            rL * H * (1 + alphaL * y(3)) - dL * y(2);
            rE * G * y(1)^2 / (KmE^2 + y(1)^2) - dE * y(3);
            rP * G * y(2)^2 / (KmP^2 + y(2)^2) - dP * y(4)];
end