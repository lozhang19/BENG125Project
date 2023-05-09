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