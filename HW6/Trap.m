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
    num_steps = (t_final - t0) / h; % Ensure an integer number of steps

    % Initialize arrays for time and solution
    t_numerical = linspace(t0, t_final, num_steps);
    u_numerical = zeros(size(t_numerical));
    du_dt_numerical = zeros(size(t_numerical));
   
    % Set initial conditions
    u_numerical(1) = u0;
    du_dt_numerical(1) = du_dt_0;

    % Trapezoidal method
    for j = 2:num_steps
        f_n = -2 * u_numerical(j - 1); % Function value at time step n
        f_np1 = -2 * (u_numerical(j - 1) + h * du_dt_numerical(j - 1)); % Function value at time step n+1
        
        u_numerical(j) = u_numerical(j - 1) + 0.5 * h * (du_dt_numerical(j - 1) + du_dt_numerical(j));
        du_dt_numerical(j) = du_dt_numerical(j - 1) + 0.5 * h * (f_n + f_np1);
    end

    % Analytical solution for the same time vector as the numerical solution
    u_analytical = cos(sqrt(2) * t_numerical);
    du_dt_exact = -sqrt(2) * sin(sqrt(2) * t_numerical);

    % Calculate the global truncation error (GTE)
    GTE_values(i) = abs(du_dt_numerical(end) - du_dt_exact(end));

    % Plot du/dt versus u
    figure;
    plot(u_numerical, du_dt_numerical, '-o');
    hold on;
    plot(u0, du_dt_0, 'ro', 'MarkerSize', 8);
    plot(u_analytical, du_dt_exact, 'k-', 'LineWidth', 1.5); % Plot exact derivative
    title(sprintf('Trapezoidal Method - Step Size: %0.3f', h));
    xlabel('u');
    ylabel('du/dt');
    legend('Solution', 'Initial Condition', 'Exact Derivative');

    % Calculate and store slope for log-log plot
    slopes(i) = log(GTE_values(i)) / log(h); 

end

% Plot GTE versus step size on a log-log scale
figure;
loglog(step_sizes, GTE_values, 'o-');
title('GTE vs Step Size at Last Solution Point');
xlabel('Step Size');
ylabel('GTE');
grid on;

% Display slopes
disp('Slopes for log-log plot:');
disp(slopes);
