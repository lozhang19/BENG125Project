function dE2Hdt = E2_equation(t, E2, FSH)
    % Parameters
    rE = 0.0018;
    dE = 2.5;
    KmE = 3.5;

    % Equation
    dE2Hdt = rE * FSH^2/(KmE^2+FSH^2) - dE * E2(t);
end