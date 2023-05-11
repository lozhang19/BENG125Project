function dLHdt = LH_equation_sensitivity_e1(t, LH, E2, P4, rL, dL, e1)
    E2_value = E2(t, e1);
    P4_value = P4(t);
    dLHdt = rL * (E2_value - P4_value) - dL * LH;
end
