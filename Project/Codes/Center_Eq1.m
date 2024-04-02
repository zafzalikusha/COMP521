clc;
clear;
close all;

%% The Main Function
f = @(x) exp(x);

% Interval and different grid sizes
a = 0;
b = 1;
grid_sizes = [0.1, 0.05, 0.025, 0.0125];

% Initialize arrays to store errors
errors_first_derivative = zeros(size(grid_sizes));
errors_second_derivative = zeros(size(grid_sizes));

num_iterations = 1;  % Number of iterations for time measurement

% Initialize arrays to store mean times
mean_times_first_derivative = zeros(size(grid_sizes));
mean_times_second_derivative = zeros(size(grid_sizes));

% Initialize cell arrays to store x and labels for each grid size
x_values = cell(size(grid_sizes));
labels_first_derivative = cell(size(grid_sizes));
labels_second_derivative = cell(size(grid_sizes));

% Initialize cell arrays to store numerical and analytical results
numerical_results_first_derivative = cell(size(grid_sizes));
numerical_results_second_derivative = cell(size(grid_sizes));
analytical_results_first_derivative = cell(size(grid_sizes));
analytical_results_second_derivative = cell(size(grid_sizes));

% Initialize arrays to store times for each iteration
times_first_derivative = zeros(1, num_iterations);
times_second_derivative = zeros(1, num_iterations);
    
%% For Loop for different grid sizes
for i = 1:length(grid_sizes)
    h = grid_sizes(i);
    x = a:h:b;
    x_values{i} = x;
    dx = 1 / h^2;
    
    f_prime_exact = exp(x);
    f_double_prime_exact = exp(x);

    for iter = 1:num_iterations    
    tic;
    % Calculate the centered finite difference approximation of the first derivative
    f_prime_approx = zeros(size(x));
    
    for j = 2:(length(x) - 1)
        f_prime_approx(j) = (f(x(j + 1)) - f(x(j - 1))) / (2 * h);
    end
    
    times_first_derivative(iter) = toc;
    % Calculate the centered finite difference approximation of the second derivative
    f_double_prime_approx = zeros(size(x));
    
    tic;
    for j = 2:(length(x) - 1)
        f_double_prime_approx(j) = (f(x(j + 1)) - 2 * f(x(j)) + f(x(j - 1))) * dx;
    end
    times_second_derivative(iter) = toc;
    end

    % Calculate the mean times
    mean_times_first_derivative(i) = mean(times_first_derivative);
    mean_times_second_derivative(i) = mean(times_second_derivative);

    % Calculate the errors
    errors_first_derivative(i) = mean(abs(f_prime_exact(2:end - 1) - f_prime_approx(2:end - 1)));
    errors_second_derivative(i) = mean(abs(f_double_prime_exact(2:end - 1) - f_double_prime_approx(2:end - 1)));
    
    % Store numerical and analytical results
    numerical_results_first_derivative{i} = f_prime_approx;
    numerical_results_second_derivative{i} = f_double_prime_approx;
    analytical_results_first_derivative{i} = f_prime_exact;
    analytical_results_second_derivative{i} = f_double_prime_exact;
    
    % Create labels with grid size and error for both derivatives
    labels_first_derivative{i} = sprintf('h = %.4f, Error (1st Derivative) = %.6f', h, errors_first_derivative(i));
    labels_second_derivative{i} = sprintf('h = %.4f, Error (2nd Derivative) = %.6f', h, errors_second_derivative(i));
end
%% Plotting Mean Times
figure;

% Plot for the mean time of the first derivative
subplot(2, 1, 1);
plot(grid_sizes, mean_times_first_derivative, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Grid Size (h)');
ylabel('Mean Time (s)');
title('Mean Time Analysis: First Derivative');
grid on;
% Add text annotations for mean time values
text(grid_sizes, mean_times_first_derivative, num2str(mean_times_first_derivative', '%.8f'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot for the mean time of the second derivative
subplot(2, 1, 2);
plot(grid_sizes, mean_times_second_derivative, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Grid Size (h)');
ylabel('Mean Time (s)');
title('Mean Time Analysis: Second Derivative');
grid on;
% Add text annotations for mean time values
text(grid_sizes, mean_times_second_derivative, num2str(mean_times_second_derivative', '%.8f'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%% Perform linear regression to find out the error order
log_grid_sizes = log(grid_sizes);
errors_first_derivative_log = log(errors_first_derivative);

% Fit a linear model
coeffs1 = polyfit(log_grid_sizes, errors_first_derivative_log, 1);

% Extract the slope (coefficient) from the model
slope1 = coeffs1(1);

disp(['The slope of the error plot for the first derivative is approximately ', num2str(slope1)]);

log_grid_sizes = log(grid_sizes);
errors_second_derivative_log = log(errors_second_derivative);

% Fit a linear model
coeffs2 = polyfit(log_grid_sizes, errors_second_derivative_log, 1);

% Extract the slope (coefficient) from the model
slope2 = coeffs2(1);

disp(['The slope of the error plot for the second derivative is approximately ', num2str(slope2)]);
%% Plot the log-log plots of errors versus grid sizes
figure;

% Plot for the first derivative
loglog(grid_sizes, errors_first_derivative, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;

% Plot for the second derivative
loglog(grid_sizes, errors_second_derivative, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');

xlabel('Grid Size (h)');
ylabel('Mean Error');
title('Error Analysis: Log-Log Plot');
legend('1st Derivative', '2nd Derivative', 'Location', 'Best');
grid on;