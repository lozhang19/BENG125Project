% Simulation parameters
startTime = 0;
stopTime = 56;
dt = 0.001;

% Parameters for E2 and Ih functions
e0 = 37.4;
e1 = 150;
e2 = 12;
e3 = 115;
e4 = 18;

% E2 function (from provided data)
E2 = @(t) e0 + e1 * exp(-(t-14)^2/e2) + e3 * exp(-(t-23)^2/e4) +e1 * exp(-(t-45)^2/e2) + e3 * exp(-(t-54)^2/e4);

% Parameters for Ih function
i0 = 0.4;
i1 = 4;
i2 = 6;

% Ih function
Ih = @(t) i0 + i1 * exp(-(t-22)^2/i2)+i2*exp(-(t-53)^2/i2);

% Initial condition for FSH
initialFSH = 14;

% Run the simulation
[t, FSH] = ode45(@(t, FSH) FSH_equation(t, FSH, E2, Ih), [startTime, stopTime], initialFSH);

% Plot results
figure;
plot(t, FSH, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('FSH Concentration');
title('1D Hormonal Regulation Model: FSH');

% Calculate E2 and Ih for each time point
E2_values = arrayfun(E2, t);
Ih_values = arrayfun(Ih, t);

% Plot E2
figure;
plot(t, E2_values, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('E2 Concentration');
title('1D Hormonal Regulation Model: E2');

%% Fixed Points
% Parameters
rF = 1.5;
E2_value = E2(14); 
Ih_value = Ih(14);  

% Use fsolve to find the fixed point
FSH_fixed_point = fsolve(@(FSH) FSH_steady_state(FSH, E2_value, Ih_value, rF), initialFSH);

disp(['FSH fixed point: ', num2str(FSH_fixed_point)]);

%% Sensitivity_gammaF

% Vary the rF parameter
gammaF_values = linspace(-3, 3, 10);  % Create a range of gammaF values to test

% Initialize a matrix to store FSH concentrations
FSH_values = zeros(length(gammaF_values), length(t));

% Loop through each rF value
for i = 1:length(gammaF_values)
    % Update the FSH_equation function with the new rF value
    FSH_equation_new = @(t, FSH) FSH_equation_sensitivity(t, FSH, E2, Ih, gammaF_values(i));
    
    % Run the simulation
    [t_temp, FSH_temp] = ode45(@(t, FSH) FSH_equation_new(t, FSH), [startTime, stopTime], initialFSH);
    
    % Interpolate FSH_temp to match the size of t
    FSH_temp_interp = interp1(t_temp, FSH_temp, t);
    
    % Store the FSH concentrations for the current rF value
    FSH_values(i, :) = FSH_temp_interp';
end

% Plot the sensitivity analysis results
figure;
surf(t, gammaF_values, FSH_values);
xlabel('Time (days)');
ylabel('gammaF Parameter');
zlabel('FSH Concentration');
title('Sensitivity Analysis: FSH vs gammaF Parameter');

% Create a new figure
figure;

% First subplot: Time vs FSH concentration
subplot(2, 1, 1);
for i = 1:length(gammaF_values)
    plot(t, FSH_values(i, :), 'LineWidth', 2);
    hold on;
end

xlabel('Time (days)');
ylabel('FSH Concentration');
title('Sensitivity Analysis: Time vs FSH Concentration');
legend(arrayfun(@(x) sprintf('gammaF = %.2f', x), gammaF_values, 'UniformOutput', false), 'Location', 'bestoutside');

% Second subplot: gammaF vs FSH concentration at specific time points
subplot(2, 1, 2);
specific_time_points = [0, 7, 14, 28]; % Set the desired time points
num_time_points = length(specific_time_points);
time_idx = zeros(1, num_time_points);

% Find the indices corresponding to the specific time points
for i = 1:num_time_points
    [~, time_idx(i)] = min(abs(t - specific_time_points(i)));
end

% Plot each time point in the same subplot with a unique color and linestyle
plot_styles = {'-', '-', '-', '-'};
colors = {'b', 'r', 'g', 'k'};
for i = 1:num_time_points
    plot(gammaF_values, FSH_values(:, time_idx(i)), 'LineWidth', 2, ...
        'Color', colors{i}, 'LineStyle', plot_styles{i});
    hold on;
end
xlabel('gammaF Parameter');
ylabel('FSH Concentration');

% Create the legend with the time points
legend_labels = cell(1, num_time_points);
for i = 1:num_time_points
    legend_labels{i} = sprintf('t = %.1f days', t(time_idx(i)));
end
legend(legend_labels, 'Location', 'best');

% Set the title
title('Sensitivity Analysis: gammaF vs FSH Concentration at Specific Time Points');

%% Sensitivity with respect to FSH

% Vary the FSH parameter
FSH_values = linspace(0, 100, 10);  % Create a range of FSH values to test

% Initialize a matrix to store FSH concentrations
FSH_sensitivity = zeros(length(FSH_values), length(t));

% Loop through each FSH value
for i = 1:length(FSH_values)
    % Run the simulation with the new initial FSH value
    [t_temp, FSH_temp] = ode45(@(t, FSH) FSH_equation(t, FSH, E2, Ih), [startTime, stopTime], FSH_values(i));
    
    % Interpolate FSH_temp to match the size of t
    FSH_temp_interp = interp1(t_temp, FSH_temp, t);
    
    % Store the FSH concentrations for the current FSH value
    FSH_sensitivity(i, :) = FSH_temp_interp';
end

% Plot the sensitivity analysis results
figure;
surf(t, FSH_values, FSH_sensitivity);
xlabel('Time (days)');
ylabel('FSH Parameter');
zlabel('FSH Concentration');
title('Sensitivity Analysis: FSH vs FSH Parameter');

% Create a new figure
figure;

% First subplot: Time vs FSH concentration
subplot(2, 1, 1);
for i = 1:length(FSH_values)
    plot(t, FSH_sensitivity(i, :), 'LineWidth', 2);
    hold on;
end

xlabel('Time (days)');
ylabel('FSH Concentration');
title('Sensitivity Analysis: Time vs FSH Concentration');
legend(arrayfun(@(x) sprintf('FSH = %.2f', x), FSH_values, 'UniformOutput', false), 'Location', 'bestoutside');

% Second subplot: FSH vs FSH concentration at specific time points
subplot(2, 1, 2);
specific_time_points = [0, 7, 14, 28]; % Set the desired time points
num_time_points = length(specific_time_points);
time_idx = zeros(1, num_time_points);

% Find the indices corresponding to the specific time points
for i = 1:num_time_points
    [~, time_idx(i)] = min(abs(t - specific_time_points(i)));
end

% Plot each time point in the same subplot with a unique color and linestyle
plot_styles = {'-', '-', '-', '-'};
colors = {'b', 'r', 'g', 'k'};
for i = 1:num_time_points
    plot(FSH_values, FSH_sensitivity(:, time_idx(i)), 'LineWidth', 2, ...
        'Color', colors{i}, 'LineStyle', plot_styles{i});
    hold on;
end
xlabel('FSH Parameter');
ylabel('FSH Concentration');

% Create the legend with the time points
legend_labels = cell(1, num_time_points);
for i = 1:num_time_points
    legend_labels{i} = sprintf('t = %.1f days', t(time_idx(i)));
end
legend(legend_labels, 'Location', 'best');

% Set the title
title('Sensitivity Analysis: FSH vs FSH Concentration at Specific Time Points');
