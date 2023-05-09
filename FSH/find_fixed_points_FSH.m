function fixedPoints = find_fixed_points_FSH(E2, Ih)
    % Define the equation to solve for the fixed point
    fixedPointEquation = @(FSH) FSH_equation(0, FSH, E2, Ih);

    % Set options for fsolve
    options = optimoptions('fsolve', 'Display', 'none');

    % Find the fixed point(s) using fsolve
    fixedPoints = fsolve(fixedPointEquation, 1, options);
end