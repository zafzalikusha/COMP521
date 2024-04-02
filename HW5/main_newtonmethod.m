clc;
clear;
close all;

% Define the functions f(x) and its derivatives (if available)
f_a = @(x) -x.^2 + x + 2;
df_a = @(x) -2.*x + 1;

f_b = @(x) exp(x) - 2 - x;
df_b = @(x) exp(x) - 1;

% Set the tolerance
TolX = 1e-9;

% Set the maximum number of iterations
MaxIter = 100;

% Define the starting points p0
p0 = [-3.0, 0.0, 6.0];

% Initialize arrays to store results
x_results_a = zeros(1, length(p0));
x_results_b = zeros(1, length(p0));

% Initialize arrays to store the number of iterations for convergence
iter_a = zeros(1, length(p0));
iter_b = zeros(1, length(p0));

% Iterate through the starting points and apply Newton's method for (a) and (b)
for i = 1:length(p0)
% For (a)
x0 = p0(i);
% Save the iteration results, which we can see how close it gets each step
[x_results_a(i), ~, function_a] = newton(f_a, df_a, x0, TolX, MaxIter);

% Plot of the results
figure
plot(function_a);
title([ 'Equation 1 starting at ', num2str(x0)]);

% Log plot of errors
true1 = 2; true2 = -1;
figure
if ( abs(function_a(end) - 2 )< 1e-8)
loglog( abs(function_a - 2) ); % true1 !!
else
loglog( abs(function_a + 1) ); % true 2
end
title(['Newton Error per step starting at ', num2str(x0)]);
xlabel('Step(a)'); 
ylabel('Error(a)');

%For (b)
x0 = p0(i);
%Save the iteration results, which we can see how close it gets each step
[x_results_b(i), ~, function_b] = newton(f_b, df_b, x0, TolX, MaxIter);

%Plot of the results
figure
plot(function_b);
title([ 'Equation 2 starting at ', num2str(x0)]);

%Log plot of errors
true1 = 2; true2 = -1;
figure
if ( abs(function_b(end) - 2 )< 1e-8)
loglog( abs(function_b - 2) );
else
loglog( abs(function_b + 1) );
end
title(['Newton Error per step starting at ', num2str(x0)]);
xlabel('Step(b)');
ylabel('Error(b)');
end

% Display the results
disp('For function f(x)_a = -x^2 + x + 2:');
disp('Starting points p0:');
disp(p0);
disp('Approximate solutions for (a):');
disp(x_results_a);

disp('For function f(x)_b = exp(x) - 2 - x:');
disp('Starting points p0:');
disp(p0);
disp('Approximate solutions for (b):');
disp(x_results_b);

% Plot function f_a
%% This needed to go farther
x_values_a = linspace(-20, 20, 1000); % Define x range for plotting function f_a
y_values_a = f_a(x_values_a); % Calculate y values for f_a
figure;
plot(x_values_a, y_values_a, 'LineWidth', 1.5); % Plot f_a
hold on;
scatter(x_results_a, f_a(x_results_a), 'ro', 'filled'); % Plot roots of f_a
title('Function f(x) = -x^2 + x + 2 and Its Roots');
xlabel('x');
ylabel('f(x)');
legend('f_a(x) = -x^2 + x + 2', 'Roots of f_a');

% Plot function f_b
x_values_b = linspace(-3, 3, 1000); % Define x range for plotting function f_b
y_values_b = f_b(x_values_b); % Calculate y values for f_b
figure;
plot(x_values_b, y_values_b, 'LineWidth', 1.5); % Plot f_b
hold on;
scatter(x_results_b, f_b(x_results_b), 'ro', 'filled'); % Plot roots of f_b
title('Function f(x) = exp(x) - 2 - x and Its Roots');
xlabel('x');
ylabel('f(x)');
legend('f_b(x) = exp(x) - 2 - x', 'Roots of f_b');

figure
errors_a = abs(function_a - (2));
loglog((errors_a(1:end-1)), (errors_a(2:end)), 'bo');
s_a = polyfit(log(errors_a(1:end-1)), log((errors_a(1:end-1))),1);
hold on
x_a = (errors_a(1:end-1));
loglog( x_a ,  abs(s_a(2)*10^12).*x_a.^s_a(1) , 'r-');
title('Convergnece Rate f(a) @ starting point = 6');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off

figure
errors_a = abs(function_a - (-1));
loglog((errors_a(1:end-1)), (errors_a(2:end)), 'bo');
s_a = polyfit(log(errors_a(1:end-1)), log((errors_a(1:end-1))),1);
hold on
x_a = (errors_a(1:end-1));
loglog( x_a ,  abs(s_a(2)*10^15).*x_a.^s_a(1) , 'r-');
title('Convergnece Rate f(a) @ starting point = -3, 0');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off

figure
errors_b = abs(function_b - (-1.841405660436961 ));
loglog((errors_b(1:end-1)), (errors_b(2:end)), 'bo');
s_b = polyfit(log(errors_b(1:end-1)), log((errors_b(1:end-1))),1);
hold on
x_b = (errors_b(1:end-1));
loglog( x_b ,  abs(s_b(2)*2*(10^15)).*x_b.^s_b(1) , 'r-');
title('Convergnece Rate f(b) @ starting point = -3');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off

figure
errors_b = abs(function_b - (1.146193220620582 ));
loglog((errors_b(1:end-1)), (errors_b(2:end)), 'bo');
s_b = polyfit(log(errors_b(1:end-1)), log((errors_b(1:end-1))),1);
hold on
x_b = (errors_b(1:end-1));
loglog( x_b ,  abs(s_b(2)*10^15).*x_b.^s_b(1) , 'r-');
title('Convergnece Rate f(b) @ starting point = 6');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off