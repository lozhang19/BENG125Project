% Parameters
rF = 0.192;
rL = 0.192;
rE = 0.16;
rP = 0.16;
dF = 0.4;
dL = 0.4;
dE = 0.4;
dP = 0.4;
alphaF = 0.2;
alphaL = 0.2;
KmE = 0.8;
KmP = 0.8;
tauH = 28;
tauG = 28;

% Hormone functions
H = @(t) 1 + 0.5 * sin(2 * pi * t / tauH);
G = @(t) 1 + 0.5 * sin(2 * pi * t / tauG);

% ODEs
dFSHdt = @(t, FSH, LH, E2, P4) rF * H(t) * (1 + alphaF * E2) - dF * FSH;
dLHdt = @(t, FSH, LH, E2, P4) rL * H(t) * (1 + alphaL * E2) - dL * LH;
dE2dt = @(t, FSH, LH, E2, P4) rE * G(t) * FSH^2 / (KmE^2 + FSH^2) - dE * E2;
dP4dt = @(t, FSH, LH, E2, P4) rP * G(t) * LH^2 / (KmP^2 + LH^2) - dP * P4;

% System of ODEs
hormone_system = @(t, y) [dFSHdt(t, y(1), y(2), y(3), y(4));
                          dLHdt(t, y(1), y(2), y(3), y(4));
                          dE2dt(t, y(1), y(2), y(3), y(4));
                          dP4dt(t, y(1), y(2), y(3), y(4))];

% Initial conditions and time span
initial_conditions = [1, 1, 1, 1];
t_span = [0, 28];

% Solve ODEs
[t, y] = ode45(hormone_system, t_span, initial_conditions);

% Plot results
figure;
subplot(2, 2, 1);
plot(t, y(:, 1));
xlabel('Time (days)');
ylabel('FSH');
title('FSH');

subplot(2, 2, 2);
plot(t, y(:, 2));
xlabel('Time (days)');
ylabel('LH');
title('LH');

subplot(2, 2, 3);
plot(t, y(:, 3));
xlabel('Time (days)');
ylabel('E2');
title('E2');

subplot(2, 2, 4);
plot(t, y(:, 4));
xlabel('Time (days)');
ylabel('P4');
title('P4');

sgtitle('Hormone Concentrations in the Female Reproductive Cycle');
