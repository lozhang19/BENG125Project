%% 2D
clear all; close all;

% Simulation parameters
startTime = 0;
stopTime = 28;
dt = 0.001;

% E2 and Ih functions parameters
e0 = 37.4;
e1 = 150;
e2 = 12;
e3 = 115;
e4 = 18;

i0 = 0.4;
i1 = 4;
i2 = 6;

E2 = @(t) e0 + e1 * exp(-(t-14).^2/e2) + e3 * exp(-(t-23).^2/e4) + e1 * exp(-(t-45).^2/e2) + e3 * exp(-(t-54).^2/e4);
Ih = @(t) i0 + i1 * exp(-(t-22)^2/i2) + i2 * exp(-(t-53)^2/i2);

% Initial conditions
Y0 = [14; 0];

% Time vector
t = startTime:dt:stopTime;

% Different gammaF values
gammaFs = [0.05, 0.1, 0.15, 0.2, 0.25];

% Prepare 3D plot
figure;

% Loop over gammaF values
for gammaF = gammaFs
    % Run RK4 method
    [Y, ~] = ODE_RK4(t, dt, Y0, E2, Ih, gammaF);
    
    % Plot in 3D
    plot3(Y(1,:), Y(2,:), t, 'LineWidth', 2);
    hold on;
end

% Set up the 3D plot
xlabel('FSH Concentration');
ylabel('E2 Concentration');
zlabel('Time (days)');
title('2D Hormonal Regulation Model: 3D Trajectory');
legend(arrayfun(@(gammaF) ['GammaF = ', num2str(gammaF)], gammaFs, 'UniformOutput', false));
hold off;

% Find fixed points on day 7, 14, 21 for gammaF = 0.1
gammaF = 0.1;
for day = [7, 14, 21]
    fixed_point = fsolve(@(Y) dydt(day, Y, E2, Ih, gammaF), Y0);
    fprintf('Fixed point for gammaF = %.2f on day %d: FSH = %.2f, E2 = %.2f\n', gammaF, day, fixed_point(1), fixed_point(2));
end

% System of equations
function [dYdt, Output] = dydt(t, Y, E2, Ih, gammaF)
    KmE = 3.5;
    dy1 = 1.5 * ( 1.3 * E2(t)) * (1 - (Ih(t)/(gammaF + Ih(t)))) - 0.4 * Y(1);
    dy2 = 0.0018 * Y(1)^2/(KmE^2+Y(1)^2) - 2.5 * Y(2);
    dYdt = [dy1; dy2];
    Output = [];
end

% RK4 method
function [Y, Output] = ODE_RK4(t, h, Y0, E2, Ih, gammaF)
    Y = zeros(length(Y0), length(t));
    Y(:,1) = Y0;
    for ii = 1:length(t)-1
        yn = Y(:,ii);
        tn = t(ii);
        K1 = dydt(tn, yn, E2, Ih, gammaF);
        K2 = dydt(tn+h/2, yn+h/2*K1, E2, Ih, gammaF);
        K3 = dydt(tn+h/2, yn+h/2*K2, E2, Ih, gammaF);
        K4 = dydt(tn+h, yn+h*K3, E2, Ih, gammaF);
        Y(:,ii+1) = yn + h/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    Output = [];
end