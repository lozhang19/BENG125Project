% Simulation parameters
startTime = 0;
stopTime = 56;
dt = 0.001;

% Parameters for E2 and Ih functions
e0 = 62.48;
e1 = 230;
e2 = 5;
e3 = 115;
e4 = 20;

% E2 function
E2 = @(t) e0 + e1 * exp(-(t-14)^2/e2)+e1 * exp(-(t-45)^2/e2);

% Parameters for E2 and Ih functions
p0= 30;
p1= 15;
p2= 100;

% P4 function
P4 = @(t) p0+ p1 * exp(-(t-22)^2/p2)+p1 * exp(-(t-53)^2/p2);

% Initial condition for FSH
initialLH = 20;

% Run the simulation
[t, LH] = ode45(@(t, LH) LH_equation(t, LH, E2, P4), [startTime, stopTime], initialLH);

% Plot results
figure;
plot(t, LH, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('LH Concentration');
title('1D Hormonal Regulation Model: LH');

% Calculate P4 for each time point
P4_values = arrayfun(P4, t);

% Plot E2
figure;
plot(t, P4_values, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('P4 Concentration');
title('1D Hormonal Regulation Model: P4');

%% Fixed Points
% Parameters
rL = 1.5;
E2_value = E2(14);  
P4_value = P4(14);  

% Use fsolve to find the fixed point
LH_fixed_point = fsolve(@(LH) LH_steady_state(LH, E2_value, P4_value, rL), initialLH);

disp(['LH fixed point: ', num2str(LH_fixed_point)]);

%% Sensitivity

% Vary the rL parameter
rL_values = linspace(0.5, 3, 10);  % Create a range of rF values to test

% Initialize a matrix to store FSH concentrations
LH_values = zeros(length(rL_values), length(t));

% Loop through each rF value
for i = 1:length(rL_values)
    % Update the FSH_equation function with the new rF value
    LH_equation_new = @(t, LH) LH_equation_sensitivity(t, LH, E2, P4, rL_values(i));
    
    % Run the simulation
    [t_temp, LH_temp] = ode45(@(t, LH) LH_equation_new(t, LH), [startTime, stopTime], initialLH);
    
    % Interpolate FSH_temp to match the size of t
    LH_temp_interp = interp1(t_temp, LH_temp, t);
    
    % Store the FSH concentrations for the current rF value
    LH_values(i, :) = LH_temp_interp';
end

% Plot the sensitivity analysis results
figure;
surf(t, rL_values, LH_values);
xlabel('Time (days)');
ylabel('rL Parameter');
zlabel('LH Concentration');
title('Sensitivity Analysis: LH vs rL Parameter');
%% Sensitivity for alphaL

% Vary the rL parameter
dL_values = linspace(-2, 2, 10);  % Create a range of rF values to test

% Initialize a matrix to store FSH concentrations
LH_values = zeros(length(rL_values), length(t));

% Loop through each rF value
for i = 1:length(rL_values)
    % Update the FSH_equation function with the new rF value
    LH_equation_new = @(t, LH) LH_equation_sensitivity(t, LH, E2, P4, dL_values(i));
    
    % Run the simulation
    [t_temp, LH_temp] = ode45(@(t, LH) LH_equation_new(t, LH), [startTime, stopTime], initialLH);
    
    % Interpolate FSH_temp to match the size of t
    LH_temp_interp = interp1(t_temp, LH_temp, t);
    
    % Store the FSH concentrations for the current rF value
    LH_values(i, :) = LH_temp_interp';
end

% Plot the sensitivity analysis results
figure;
surf(t, dL_values, LH_values);
xlabel('Time (days)');
ylabel('alphaL Parameter');
zlabel('LH Concentration');
title('Sensitivity Analysis: LH vs dL Parameter');