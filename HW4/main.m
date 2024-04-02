% Homework 04
% Use numerical integration to approximate the integral of a function f(x)
% over the interval x \in [a, b]
close all;
clear;
clc;

% Define the function
Fx = @(x) exp(-2 .* x) .* sin(2 .* pi .* x);

a = 0;
b = 3.5;

% Calculate the exact integral
exact_integral = (exp(-2 * b) - exp(-2 * a)) / 2 + (sin(2 * pi * a) - sin(2 * pi * b)) / (4 * pi^2);

% Define the exact integral value
exact_integral = 0.1446445197;

% Define the number of subintervals
N_values = [20, 40, 80, 160];

% Initialize arrays to store errors, orders of accuracy, and iterations
errors_trapezoidal = zeros(size(N_values));
orders_trapezoidal = zeros(size(N_values));
iterations_trapezoidal = zeros(size(N_values));

errors_simpson = zeros(size(N_values));
orders_simpson = zeros(size(N_values));
iterations_simpson = zeros(size(N_values));

%% Loop through different values of N
for k = 1:length(N_values)
    N = N_values(k);
    
    % Part 1
    % Apply the composite trapezoidal rule to calculate the integral
    a = 0;
    b = 3.5;
    
    % Initialize a counter for iterations
    counter_trapezoidal = 0;
    
    % Perform the trapezoidal integration
    for i = 1:N
        % Calculate Itc using the trapezoidal rule
        [Itc] = traprl(Fx, a, b, i);
        
        % Count the iterations
        counter_trapezoidal = counter_trapezoidal + 1;
        
        % Calculate the absolute error
        error_trapezoidal = abs(exact_integral - Itc);
        errors_trapezoidal(k) = error_trapezoidal;
        
        % Calculate the order of accuracy
        if k > 1
            orders_trapezoidal(k) = -diff(log(errors_trapezoidal(k-1:k))) / log(2);
        end
    end
    
    % Store the number of iterations for the trapezoidal method
    iterations_trapezoidal(k) = counter_trapezoidal;
    
    %% Part 2
    % Apply the composite Simpson's Rule to calculate the integral
    % Initialize a counter for iterations
    counter_simpson = 0;
    
    % Perform the Simpson's integration
    for i = 1:N % Starting from 2 to utilize the composite Simpson's Rule
        [Isc] = simprl(Fx, a, b, i);
        
        % Count the iterations
        counter_simpson = counter_simpson + 1;
        
        % Calculate the absolute error
        error_simpson = abs(exact_integral - Isc);
        errors_simpson(k) = error_simpson;
        
        % Calculate the order of accuracy
        if k > 1
            orders_simpson(k) = -diff(log(errors_simpson(k-1:k))) / log(2);
        end
    end
    
    % Store the number of iterations for the Simpson's method
    iterations_simpson(k) = counter_simpson;
end
%%
a = 0;
b = 3.5;
% Calculate the error for the trapezoidal method with N = 160
N_trap = 160;
desired_error_Trapezoidal = abs(exact_integral - traprl(Fx, a, b, N_trap));

% Initialize variables to track iterations and results for Trapezoidal Rule
iterations_adaptive_simpson_T = 0;
Ias_T = 0;

% Loop to adaptively reduce error for Trapezoidal Rule
while N_trap > 1
    [SRmat_T, Ias_T, ~] = adapt(Fx, a, b, desired_error_Trapezoidal); % Use the exact error
    iterations_adaptive_simpson_T = size(SRmat_T, 1);
    
    % Check if the error is less than the desired error
    if abs(Ias_T - exact_integral) < desired_error_Trapezoidal
        break;
    end
    
    % Reduce N_trap for the next iteration
    N_trap = N_trap / 2;
end

% Calculate the error for Simpson's Rule at N = 160
N_Simpson = 160;
desired_error_Simpson = abs(exact_integral - simprl(Fx, a, b, N_Simpson));

% Initialize variables to track iterations and results for Simpson's Rule
iterations_adaptive_simpson_S = 0;
Ias_S = 0;

% Loop to adaptively reduce error for Simpson's Rule
while N_Simpson > 1
    [SRmat_S, Ias_S, ~] = adapt(Fx, a, b, desired_error_Simpson); % Use the exact error
    iterations_adaptive_simpson_S = size(SRmat_S, 1);
    
    % Check if the error is less than the desired error
    if abs(Ias_S - exact_integral) < desired_error_Simpson
        break;
    end
    
    % Reduce N_Simpson for the next iteration
    N_Simpson = N_Simpson / 2;
end
%% Display the number of iterations needed for Simpson's Rule
fprintf('Number of iterations (Adaptive Simpson''s Rule) for error %.16f (Simpson''s Rule N=%d): %d\n', desired_error_Simpson, N_Simpson, iterations_adaptive_simpson_S);

% Display the number of iterations needed for Trapezoidal Rule
fprintf('Number of iterations (Adaptive Simpson''s Rule) for error %.16f (Trapezoidal Rule N=%d): %d\n', desired_error_Trapezoidal, N_trap, iterations_adaptive_simpson_T);

% Display errors and orders of accuracy
fprintf('Errors (Trapezoidal Rule):\n');
fprintf('%.15f\n', errors_trapezoidal);

disp('Orders of Accuracy (Trapezoidal Rule):');
disp(orders_trapezoidal);

fprintf('Errors (Simpson Rule):\n');
fprintf('%.15f\n', errors_simpson);

disp('Orders of Accuracy (Simpson''s Rule):');
disp(orders_simpson);

%% Define the data for creating the table
N = [20; 40; 80; 160];
trapezoidal_error = [0.016275064552135; 0.004026532514390; 0.001003994100437; 0.000250833652897];
trapezoidal_order = [NaN; 2.0151; 2.0038; 2.0009];
simpson_error = [0.000056311498191; 0.000003518704214; 0.000000219829617; 0.000000013713331];
simpson_order = [NaN; 4.0003; 4.0006; 4.0027];

% Create a table
T = table(N, trapezoidal_error, trapezoidal_order, simpson_error, simpson_order);

% Define the variable names (table headers)
T.Properties.VariableNames = {'N', 'Trapezoidal_Error', 'Trapezoidal_Order', 'Simpson_Error', 'Simpson_Order'};

% Display the table
disp(T);

% Plot the errors for Trapezoidal and Simpson's Rule
figure;
semilogy(N_values, errors_trapezoidal, '-o', 'DisplayName', 'Trapezoidal Rule');
hold on;
semilogy(N_values, errors_simpson, '-o', 'DisplayName', "Simpson's Rule");

% Add labels and legend
xlabel('Number of Subintervals (N)');
ylabel('Error (log scale)');
title('Convergence of Trapezoidal and Simpson''s Rule');
legend('Location', 'Best');
grid on;

% Customize the plot
set(gca, 'XTick', N_values);
set(gca, 'XTickLabel', arrayfun(@num2str, N_values, 'UniformOutput', false));
hold off;