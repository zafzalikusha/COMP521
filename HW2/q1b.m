clc;
clear;
close all;
%% The Main Function
f = @(x) 2 * (cos(pi * x)).^2 - 1;

%  The interval and different grid sizes
a = -1;
b = 1;
grid_sizes = [0.02, 0.01, 0.005, 0.002];
h = grid_sizes;

% Initialize arrays to store error 
errors = zeros(size(grid_sizes));

% Create cell arrays to store x and labels for each grid size
x_values = cell(size(grid_sizes));
labels = cell(size(grid_sizes));
 
% Initialize arrays to store numerical and analytical results
numerical_results = cell(size(grid_sizes));
analytical_results = cell(size(grid_sizes));
%% For Loop for different grid sizes
for i = 1:length(grid_sizes)
    h = grid_sizes(i);
    x = a:h:b;
    x_values{i} = x;
    dx = 1/h^2;
    
    % Calculate the analytical second derivative
     f_double_prime_exact = 4 * pi^2 * ((sin(pi * x)).^2 - (cos(pi * x)).^2);
    
   % Calculate the centered finite difference approximation of the second derivative
    f_double_prime_approx = zeros(size(x));
    
    for j = 2:(length(x) - 1)
        f_double_prime_approx(j) = (f(x(j + 1)) - 2 * f(x(j)) + f(x(j - 1))) * dx;
    end 
   
    % Calculate the error
    errors(i) = max(abs(f_double_prime_exact(2:length(grid_sizes)-1) - f_double_prime_approx(2:length(grid_sizes)-1)));
    
    % Store numerical and analytical results
    numerical_results{i} = f_double_prime_approx;
    analytical_results{i} = f_double_prime_exact;
    
    % Create labels with grid size and error
    labels{i} = sprintf('h = %.4f, Error = %.6f', h, errors(i));
end
%% Plot the log-log plots of errors versus grid sizes with labels

figure;
loglog(grid_sizes, errors, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
text(grid_sizes, errors, labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
xlabel('Grid Size (h)');
ylabel('Error');
title('Error Analysis: Log-Log Plot');
grid on;

% Plot the numerical and analytical results for h = 0.02 in the same plot
h_index = find(grid_sizes == 0.02); % Find the index for h = 0.02
figure;
plot(x_values{h_index}, numerical_results{h_index}, 'k--', 'LineWidth', 2, 'DisplayName', 'Numerical (h = 0.02)');
hold on;
plot(x_values{h_index}, analytical_results{h_index}, 'yo', 'LineWidth', 2, 'DisplayName', 'Analytical (h = 0.02)');
xlabel('x');
ylabel("f''(x)");
title('Numerical vs. Analytical Results (h = 0.02)');
legend('Location', 'northwest');
grid on;
%% Perform linear regression
log_grid_sizes = log(grid_sizes);
log_errors = log(errors);

% Fit a linear model
coeffs = polyfit(log_grid_sizes, log_errors, 1);

% Extract the slope (coefficient) from the model
slope = coeffs(1);

disp(['The slope of the error plot is approximately ', num2str(slope)]);
