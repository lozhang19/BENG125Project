function matcont_analysis(E2, Ih)
    % Define the system
    sys = @(t, FSH) FSH_equation(t, FSH, E2, Ih);

    % Set the initial point and parameter values
    initial_point = [4.5];
    parameters = [];

    % Continuation options
    opt = contset;
    opt.Singularities = 1;
    opt.Backward = 0;

    % Perform the continuation and bifurcation analysis
    [fixed_points, ~, ~, ~] = cont(@equilibrium, initial_point, parameters, opt);

    % Plot the results (example: FSH values)
    figure;
    plot(fixed_points(1, :), 'b', 'LineWidth', 2);
    xlabel('Continuation Steps');
    ylabel('FSH');
    title('Continuation and Bifurcation Analysis');
end
