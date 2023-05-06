function J = hormone_system_jacobian(t, y, rF, rL, rE, rP, dF, dL, dE, dP, alphaF, alphaL, KmE, KmP, tauH, tauG)
    H = 1 + 0.5 * sin(2 * pi * t / tauH);
    G = 1 + 0.5 * sin(2 * pi * t / tauG);
    
    J = [-dF, 0, rF * H * alphaF, 0;
         0, -dL, rL * H * alphaL, 0;
         rE * G * 2 * y(1) / (KmE^2 + y(1)^2) - rE * G * y(1)^4 / (KmE^2 + y(1)^2)^2, 0, -dE, 0;
         0, rP * G * 2 * y(2) / (KmP^2 + y(2)^2) - rP * G * y(2)^4 / (KmP^2 + y(2)^2)^2, 0, -dP];
end