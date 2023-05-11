function dFSHdt = FSH_derivatives(FSH_range, t_range, E2, Ih, dF_single)
    dFSHdt = zeros(length(FSH_range), length(t_range));
    
    for j = 1:length(FSH_range)
        for k = 1:length(t_range)
            dFSHdt(j, k) = FSH_equation_sensitivity_dF(t_range(k), FSH_range(j), E2, Ih, dF_single);
        end
    end
end