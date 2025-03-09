% Parameters
R = 1;          % Length of the pendulum (m)
g = 9.81;       % Acceleration due to gravity (m/s^2)
alpha = deg2rad(30); % Initial angle (converted to radians)
h = 0.006;     % Time step (s)
T_total = 16;  % Total simulation time (s)
n_steps = floor(T_total / h); % Number of time steps (ensure integer)

% Initialize arrays
theta = zeros(1, n_steps);
omega = zeros(1, n_steps);
time = linspace(0, T_total, n_steps);

% Initial conditions
theta(1) = alpha;
omega(1) = 0;

% Forward Euler Method
for n = 1:n_steps-1
    omega(n+1) = omega(n) - (g / R) * sin(theta(n)) * h;
    theta(n+1) = theta(n) + omega(n) * h;
end

% Plotting the results
figure;
plot(time, theta, 'b', 'LineWidth', 1.5);
hold on;

% Ideal sinusoidal solution for comparison
omega_0 = sqrt(g / R);
ideal_theta = alpha * cos(omega_0 * time);
plot(time, ideal_theta, 'r--', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Angle (rad)');
title('Pendulum Motion using Forward Euler Method');
legend('Forward Euler Computed', 'Ideal');
grid on;