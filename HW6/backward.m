% Backward Euler Method for the given ODE
clc;
close all;
clear;

% Define parameters
t0 = 0;
t_final = 10;
u0 = 1;
du_dt_0 = 0;

% Define step sizes
step_sizes = [0.1, 0.05, 0.01, 0.005];

% Initialize GTE_values to store errors
GTE_values = zeros(size(step_sizes));

% Initialize array to store slopes
slopes = zeros(size(step_sizes));

% Loop over different step sizes
for i = 1:length(step_sizes)
    h = step_sizes(i);
    num_steps = (t_final - t0) / h;

    % Run Backward Euler simulation
    t_numerical = linspace(t0, t_final, num_steps); % Adjusted to num_steps + 1
    u_numerical = zeros(size(t_numerical));
    du_dt_numerical = zeros(size(t_numerical));

    % Set initial conditions
    u_numerical(1) = u0;
    du_dt_numerical(1) = du_dt_0;

    % Backward Euler method
    for j = 2:num_steps % Adjusted loop bounds
        u_numerical(j) = (u_numerical(j - 1) + h * du_dt_numerical(j - 1)) / (1 + 2 * h);
        du_dt_numerical(j) = (u_numerical(j) - u_numerical(j - 1)) / h;
    end

    % Analytical solution for the corresponding time points
    t_analytical = linspace(t0, t_final, num_steps); % Adjusted to num_steps + 1
    u_analytical = cos(sqrt(2) * t_analytical);
    du_dt_exact = -sqrt(2) * sin(sqrt(2) * t_analytical);

    % Calculate Global Truncation Error (GTE) at the last solution point
    GTE_values(i) = abs(u_numerical(end) - u_analytical(end));

    % Plot du/dt versus u
    figure;
    plot(u_numerical, du_dt_numerical, '-o');
    hold on;
    plot(u0, du_dt_0, 'ro', 'MarkerSize', 8);
    plot(u_analytical, du_dt_exact, 'k-', 'LineWidth', 1.5); % Plot exact derivative
    title(sprintf('Backward Euler Method: Step Size = %0.3f', h));
    xlabel('u');
    ylabel('du/dt');
    legend('Numerical Solution', 'Initial Condition', 'Exact Derivative');
    grid on;
    hold off;

    % Calculate and store slope for log-log plot
    slopes(i) = log(GTE_values(i)) / log(h);
end

% Plot log-log plot of GTE versus step size
figure;
loglog(step_sizes, GTE_values, '-o');
title('Global Truncation Error vs Step Size');
xlabel('Step Size (dt)');
ylabel('Global Truncation Error (GTE)');
grid on;

% Display slopes
disp('Slopes for log-log plot:');
disp(slopes);
