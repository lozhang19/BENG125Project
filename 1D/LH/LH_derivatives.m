function dLHdt = LH_derivatives(LH_range, t_range, E2, P4)
    dLHdt = zeros(length(LH_range), length(t_range));
    
    for j = 1:length(LH_range)
        for k = 1:length(t_range)
            dLHdt(j, k) = LH_derivative(LH_range(j), t_range(k), E2, P4);
        end
    end
end
