clc;
clear;
close all;

% Function to define the differential equation dy/dx = f(x, y)
f = @(x, y) (exp(-5.*x).*sin(75.*x) + 0.0001 * rand(size(x)));

% Initial conditions
x0 = 0;
y0 = exp(-5.*x0).*sin(75.*x0);

% Time span
tspan = [0, 1];  % Adjust the end point as needed

% Different step sizes (h values)
h_values = [0.01, 0.001, 0.0001];
num_iterations = 10;

% Initialize cell arrays to store errors
errors_first_derivative = cell(size(h_values));
errors_second_derivative = cell(size(h_values));

% Initialize arrays to store max errors
mean_errors_first_derivative = zeros(size(h_values));
mean_errors_second_derivative = zeros(size(h_values));

% Initialize arrays to store mean times
mean_times_first_derivative = zeros(size(h_values));
mean_times_second_derivative = zeros(size(h_values));

% Initialize arrays to store individual iteration times
times_first_derivative = zeros(1, num_iterations);
times_second_derivative = zeros(1, num_iterations);

% Loop through different step sizes
for i = 1:length(h_values)
    h = h_values(i);
   
    % Solve the ODE using ode45 with options
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [x, y] = ode45(f, tspan, y0, options);
    
    % Calculate exact derivatives
    f_prime_exact = exp(-5.*x).*(5.*sin(75.*x)-75*cos(75.*x));
    f_double_prime_exact = exp(x);
    
    % Perform multiple iterations
    for iter = 1:num_iterations
       
        % Timing for Numerical Derivatives
        tic;
        % Calculate numerical derivatives
        y_prime_numerical = gradient(y, x);
        times_first_derivative(iter) = toc;

        tic;
        y_double_prime_numerical = gradient(y_prime_numerical, x);
        times_second_derivative(iter) = toc;
    end

 % Calculate errors including h
  errors_first_derivative{i} = h * abs(y_prime_numerical - f_prime_exact);
  errors_second_derivative{i} = h * abs(y_double_prime_numerical - f_double_prime_exact);

 % Store the mean errors
  mean_errors_first_derivative(i) = (mean(errors_first_derivative{i}));
  mean_errors_second_derivative(i) = (mean(errors_second_derivative{i}));

 % Store the mean times for both derivatives
  mean_times_first_derivative(i) = mean(times_first_derivative);
  mean_times_second_derivative(i) = mean(times_second_derivative);
end

% Plot the log-log results
figure;

% Plot for the first derivative
loglog(h_values, mean_errors_first_derivative, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'First Derivative');
hold on;

% Plot for the second derivative
loglog(h_values, mean_errors_second_derivative, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Second Derivative');

xlabel('Step Size (h)', 'FontSize', 12);
ylabel('Mean Error', 'FontSize', 12);
title('Log-Log Plot of Step Size vs. Mean Error (Including h)', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'Fontsize', 10);

% Plot mean times
figure;
subplot(2, 1, 1);
plot(h_values, mean_times_first_derivative, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Mean Time for First Derivative');
xlabel('Step Size (h)', 'FontSize', 12);
ylabel('Mean Time (s)', 'FontSize', 12);
title('Mean Time Analysis for First Derivative', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'Fontsize', 10);

subplot(2, 1, 2);
plot(h_values, mean_times_second_derivative, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Mean Time for Second Derivative');
xlabel('Step Size (h)', 'FontSize', 12);
ylabel('Mean Time (s)', 'FontSize', 12);
title('Mean Time Analysis for Second Derivative', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'Fontsize', 10);