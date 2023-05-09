% Vary the rF parameter
rF_values = linspace(0.5, 3, 10);  % Create a range of rF values to test

% Initialize a matrix to store FSH concentrations
FSH_values = zeros(length(rF_values), length(t));

% Loop through each rF value
for i = 1:length(rF_values)
    % Update the FSH_equation function with the new rF value
    FSH_equation_new = @(t, FSH) FSH_equation(t, FSH, E2, Ih, rF_values(i));
    
    % Run the simulation
    [t_temp, FSH_temp] = ode45(@(t, FSH) FSH_equation_new(t, FSH), [startTime, stopTime], initialFSH);
    
    % Interpolate FSH_temp to match the size of t
    FSH_temp_interp = interp1(t_temp, FSH_temp, t);
    
    % Store the FSH concentrations for the current rF value
    FSH_values(i, :) = FSH_temp_interp';
end

% Plot the sensitivity analysis results
figure;
surf(t, rF_values, FSH_values);
xlabel('Time (days)');
ylabel('rF Parameter');
zlabel('FSH Concentration');
title('Sensitivity Analysis: FSH vs rF Parameter');