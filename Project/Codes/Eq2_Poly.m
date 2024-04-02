close all; 
clear all; 
clc;

% Define the original function and its exact first and second derivatives
f = @(x) exp(-5.*x).*sin(75.*x);
df_exact = @(x) -exp(-5.*x).*(5.*sin(75.*x)-75*cos(75.*x));
ddf_exact = @(x) -exp(-5.*x).*(5600*sin(75.*x)+750*cos(75.*x));

% Parameters
dx_values = [0.01, 0.001, 0.0001, 0.00001];
n = 10; % Maximum degree 
num_iterations = 10;

% Initialize arrays to store mean times
mean_times_first_derivative = zeros(size(dx_values));
mean_times_second_derivative = zeros(size(dx_values));

% Initialize arrays to store errors
errors_first_derivative = zeros(size(dx_values));
errors_second_derivative = zeros(size(dx_values));

% Loop over different step sizes
for idx = 1:length(dx_values)
    dx = dx_values(idx);
    xvec = [0:dx:1];

    % Generate noise
    noise = 0.0001 * rand(size(xvec));

    % Add noise to the function
    f_noisy = f(xvec) + noise;

    % Initialize arrays to store individual iteration times
    times_first_derivative = zeros(1, num_iterations);
    times_second_derivative = zeros(1, num_iterations);

    % Loop over iterations
    for iter = 1:num_iterations
        % Timing for First Derivative
        tic;
        first_derivative = polyder(polyfit(xvec, f_noisy, n));
        time_first_derivative(iter) = toc;

        % Timing for Second Derivative
        tic;
        second_derivative = polyder(first_derivative);
        time_second_derivative(iter) = toc;
    end

    % Evaluate the derivative polynomials
    aa_first_der_vals = polyval(first_derivative, xvec);
    aa_second_der_vals = polyval(second_derivative, xvec);

    % Calculate errors
    errors_first_derivative(idx) = max(abs(df_exact(xvec) - aa_first_der_vals));
    errors_second_derivative(idx) = max(abs(ddf_exact(xvec) - aa_second_der_vals));

    % Calculate the mean times
    mean_times_first_derivative(idx) = mean(time_first_derivative);
    mean_times_second_derivative(idx) = mean(time_second_derivative);
end

% Plot the log-log plots of errors versus grid sizes
figure;

% Plot for the first derivative
loglog(dx_values, errors_first_derivative, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;

% Plot for the second derivative
loglog(dx_values, errors_second_derivative, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Grid Size (dx)');
ylabel('Max Error');
title('Error Analysis: Log-Log Plot');
legend('1st Derivative', '2nd Derivative', 'Location', 'Best');
grid on;

% Plotting Mean Times
figure;

% Plot for the mean time of the first derivative
subplot(2, 1, 1);
plot(dx_values, mean_times_first_derivative, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Step Size (dx)');
ylabel('Mean Time (s)');
title('Mean Time Analysis: First Derivative');
grid on;

% Plot for the mean time of the second derivative
subplot(2, 1, 2);
plot(dx_values, mean_times_second_derivative, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Step Size (dx)');
ylabel('Mean Time (s)');
title('Mean Time Analysis: Second Derivative');
grid on;
