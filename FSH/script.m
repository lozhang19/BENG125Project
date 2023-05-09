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