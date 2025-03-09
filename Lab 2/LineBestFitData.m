% Parameters
R = 1;          % Length of the pendulum (m)
g = 9.81;       % Acceleration due to gravity (m/s^2)
alpha = deg2rad(30); % Initial angle (converted to radians)
h = 0.006;     % Time step (s)
T_total = 16;  % Total simulation time (s)
n_steps = floor(T_total / h); % Number of time steps (ensure integer)

% Load Data
run('Data.m'); % Load Part1Data1 from Data.m

% Initialize arrays
theta = zeros(1, n_steps);
omega = zeros(1, n_steps);
time = linspace(0, T_total, n_steps);

% Initial conditions
theta(1) = alpha;
omega(1) = 0;

% 4th-Order Runge-Kutta Method
for n = 1:n_steps-1
    k1_theta = omega(n);
    k1_omega = -(g / R) * sin(theta(n));

    k2_theta = omega(n) + 0.5 * h * k1_omega;
    k2_omega = -(g / R) * sin(theta(n) + 0.5 * h * k1_theta);

    k3_theta = omega(n) + 0.5 * h * k2_omega;
    k3_omega = -(g / R) * sin(theta(n) + 0.5 * h * k2_theta);

    k4_theta = omega(n) + h * k3_omega;
    k4_omega = -(g / R) * sin(theta(n) + h * k3_theta);

    theta(n+1) = theta(n) + (h / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta);
    omega(n+1) = omega(n) + (h / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega);
end

% Elliptic Integral Calculation
k = sin(alpha / 2);  % Modulus
elliptic_integral = integral(@(theta) 1 ./ sqrt(1 - k^2 * sin(theta).^2), 0, pi/2);
T_elliptic = 4 * sqrt(R / g) * elliptic_integral;

% Line of Best Fit
exp_time = linspace(0, T_total, length(Part1Data1(1, :)));
exp_data = Part1Data1(1, :);
coeffs = polyfit(exp_time, exp_data, 1); % Linear fit [slope, intercept]
fit_line = polyval(coeffs, exp_time);

% Plotting the results
figure;
plot(time, theta, 'b', 'LineWidth', 1.5);
hold on;

% Ideal sinusoidal solution for comparison
omega_0 = sqrt(g / R);
ideal_theta = alpha * cos(omega_0 * time);
plot(time, ideal_theta, 'r--', 'LineWidth', 1.5);

% Plotting experimental data
plot(exp_time, exp_data, 'g.', 'LineWidth', 1.5);

% Plotting line of best fit
plot(exp_time, fit_line, 'k--', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Angle (rad)');
title(['Pendulum Motion using 4th-Order Runge-Kutta Method (T_{elliptic} = ' ...
    , num2str(T_elliptic, '%.4f'), ' s)']);
legend('Runge-Kutta Computed', 'Ideal', 'Experimental Data', 'Line of Best Fit');
grid on;

% Display equation of the line of best fit
equation_str = ['y = ', num2str(coeffs(1), '%.4f'), 'x + ', 
    num2str(coeffs(2), '%.4f')];
disp(['Line of Best Fit: ', equation_str]);

% Output the equation on the plot
text(0.05 * T_total, max(exp_data) * 0.9, equation_str, 'FontSize', 12, 'Color', ...
    ...
'k', 'BackgroundColor', 'w');