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

%% Sensitivity

% Vary the rF parameter
rF_values = linspace(0.5, 3, 10);  % Create a range of rF values to test

% Initialize a matrix to store FSH concentrations
FSH_values = zeros(length(rF_values), length(t));

% Loop through each rF value
for i = 1:length(rF_values)
    % Update the FSH_equation function with the new rF value
    FSH_equation_new = @(t, FSH) FSH_equation_sensitivity(t, FSH, E2, Ih, rF_values(i));
    
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
ylabel('rF Parameter');
zlabel('FSH Concentration');
title('Sensitivity Analysis: FSH vs gammaF Parameter');

%% Bifurcation
% Load the FSH_system function

system_name = 'FSH_system';

% Set the parameters for E2 and Ih functions (use the values from your original script)
E2_value = E2(14);
Ih_value = Ih(14);

% Define the parameter range for rF
rF_range = [0.5, 10];

% Initial condition for FSH
initialFSH = 14;

% Set up the bifurcation analysis
ap = 1; % index of the bifurcation parameter in the parameter vector
p = [rF_range(1), E2_value, Ih_value];
[x0, v0] = init_EP_EP(system_name, initialFSH, p, ap);

% Set the options for the continuation
opt = contset;
opt.MaxNumPoints = 100; % number of points in the continuation
opt.Singularities = 1; % detect singularities

% Perform the bifurcation analysis
[x, v, s, h, f] = cont(@equilibrium, x0, v0, opt);

% Plot the bifurcation diagram
figure;
plot(x(2,:), x(1,:), 'b');
xlabel('Parameter (rF)');
ylabel('FSH Steady-State Concentration');
title('Bifurcation Diagram');