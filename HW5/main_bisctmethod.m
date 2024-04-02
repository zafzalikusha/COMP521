clc;
clear;
close all;

% Set the intervals for (a) and (b)
interval_a = [-4.0, 1.0];
interval_b = [0.5, 3.0];

% Set the tolerance
TolX = 1e-9;

% Set the maximum number of iterations
MaxIter = 100;

% (a) Define the function f_a(x) as an anonymous function
f_a = @(x) -x.^2 + x + 2;

% (a) Call the bisection method for interval_a with f_a
[x_a, err_a, guesses_a] = bisct(f_a, interval_a(1), interval_a(2), TolX, MaxIter);
iterations_a = length(guesses_a);

% Plot the guesses vs iteration for function b
figure;
plot(1:iterations_a, guesses_a, 'o-');
title('Guesses vs Iteration for f_a');
xlabel('Iteration');
ylabel('Guesses');

figure
errors_a = abs(guesses_a - (-1.000000000116415));
loglog((errors_a(1:end-1)), (errors_a(2:end)), 'bo');
s_b = polyfit(log(errors_a(1:end-1)), log((errors_a(1:end-1))),1);
hold on
x_a = (errors_a(1:end-1));
loglog( x_a ,  abs(s_b(2)*10^15).*x_a.^s_b(1) , 'r-');
title('Convergnece Rate f(a)');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off

% (b) Define the function f_b(x) as an anonymous function
f_b = @(x) exp(x) - 2 - x;

% (b) Call the bisection method for interval_b with f_b
[x_b, err_b, guesses_b] = bisct(f_b, interval_b(1), interval_b(2), TolX, MaxIter);
iterations_b = length(guesses_b);

% Plot the guesses vs iteration for function b
figure;
plot(1:iterations_b, guesses_b, 'o-');
title('Guesses vs Iteration for f_b');
xlabel('Iteration');
ylabel('Guesses');

% Display the results for both functions
fprintf('For interval a [%f, %f] and function f(x) = -x^2 + x + 2:\n', interval_a(1), interval_a(2));
%fprintf('Approximate solution x_a: %f\n', x_a);
fprintf('Half interval width (err_a): %f\n', err_a);

fprintf('For interval b [%f, %f] and function f(x) = exp(x) - 2 - x:\n', interval_b(1), interval_b(2));
fprintf('Approximate solution x_b: %f\n', x_b);
fprintf('Half interval width (err_b): %f\n', err_b);

% Define the functions
f_a = @(x) -x.^2 + x + 2;
f_b = @(x) exp(x) - 2 - x;

% Calculate roots using the bisection method
[x_a, ~, ~] = bisct(f_a, -4.0, 1.0, 1e-9, 100);
[x_b, ~, ~] = bisct(f_b, 0.5, 3.0, 1e-9, 100);

% Define the x-values for the plots
x_values_a = linspace(-6, 6, 1000);
x_values_b = linspace(-1, 4, 1000);

% Evaluate the functions
y_values_a = f_a(x_values_a);
y_values_b = f_b(x_values_b);

% Plot f_a
figure;
plot(x_values_a, y_values_a, 'LineWidth', 1.5);
hold on;
plot(x_a, f_a(x_a), 'ro', 'MarkerSize', 8);  % Mark the root of f_a
title('Function f(x) = -x^2 + x + 2 with Root');
xlabel('x');
ylabel('f(x)');
grid on;
legend('f_a(x) = -x^2 + x + 2', 'Root of f_a');

% Plot f_b
figure;
plot(x_values_b, y_values_b, 'LineWidth', 1.5);
hold on;
plot(x_b, f_b(x_b), 'ro', 'MarkerSize', 8);  % Mark the root of f_b
title('Function f(x) = exp(x) - 2 - x with Root');
xlabel('x');
ylabel('f(x)');
grid on;
legend('f_b(x) = exp(x) - 2 - x', 'Root of f_b');

% Plot the convergence behavior for function a
figure;
semilogy(1:length(guesses_a), abs(guesses_a - x_a), 'o-');
title('Convergence Behavior of Bisection for f_a');
xlabel('Iteration');
ylabel('Absolute Error');

% Plot the convergence behavior for function b
figure;
semilogy(1:length(guesses_b), abs(guesses_b - x_b), 'o-');
title('Convergence Behavior of Bisection for f_b');
xlabel('Iteration');
ylabel('Absolute Error');

figure
errors_b = abs(guesses_b - (1.146193220163696));
loglog((errors_b(1:end-1)), (errors_b(2:end)), 'bo');
s_b = polyfit(log(errors_b(1:end-1)), log((errors_b(1:end-1))),1);
hold on
x_b = (errors_b(1:end-1));
loglog( x_b ,  abs(s_b(2)*10^15).*x_b.^s_b(1) , 'r-');
title('Convergnece Rate f(b)');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off
