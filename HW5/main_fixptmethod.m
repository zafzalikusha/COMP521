clc;
clear;
close all;

% Set the initial guesses for g_a and g_b
p0_a = [2.5; 0.15; 1.5];
p0_b = [2.5; 0.15; 0.25];

% Define function f(x), a and b
f_a = @(x) -x.^2 + x + 2;
f_b = @(x) exp(x) - 2 - x;

% Call the fixpt function for each initial guess for g_a
results_a = zeros(length(p0_a), 1);
final_errors_a = zeros(length(p0_a), 1);

for i = 1:length(p0_a)
% Give the max iterations to fixpt as the fourth argument, and save the third output,
% it called 'guessesA', we can see what the fixpt is guessing and whats happening.
[x_a, err_a, guessesA] = fixpt(f_a, p0_a(i), 1e-6, 10000);
results_a(i) = x_a;
final_errors_a(i) = err_a;
% The algorithm Progress
% figure
% plot(guessesA);
% title(['Guesses for fixpt A starting at ', num2str(p0_a(i)), ', why does this look like a brick?']);
% xlabel('Iteration');
% ylabel('Guess');
end

% Call the fixpt function for each initial guess for g_b
results_b = zeros(length(p0_b), 1);
final_errors_b = zeros(length(p0_b), 1);

for i = 1:length(p0_b)
[x_b, err_b, guessesB] = fixpt(f_b, p0_b(i), 1e-9, 10000);
results_b(i) = x_b;
final_errors_b(i) = err_b;
% figure
% plot(guessesB);
% title(['Guesses for fixpt B starting at ', num2str(p0_b(i))]);
% xlabel('Iteration');
% ylabel('Guess');
end

% Display the results for g_a
disp('Results for g_a:');
disp('Approximate fixed points:');
disp(results_a);
disp('Final errors:');
disp(final_errors_a);

% Display the results for g_b
disp('Results for g_b:');
disp('Approximate fixed points:');
disp(results_b);
disp('Final errors:');
disp(final_errors_b);

% Define the range for x values
x_values_a = linspace(-4, 4, 1000);
x_values_b = linspace(-3, 2, 1000);

% Calculate the function values for g_a and g_b
y_values_a = g_a(x_values_a);
y_values_b = g_b(x_values_b);

% Plot for g_a
% figure;
% plot(x_values_a, y_values_a, 'LineWidth', 1.5);
% xlabel('x');
% ylabel('g_a(x)');
% title('Function g_a(x) = -x^2 + x + 2');

%% Section 2
% Apply fixpt to find the root for g_a
[x_a, x_e, guessesA] = fixpt(f_a, 1.5, 1e-6, 1000000); % Change 0.5 to the desired initial guess

figure
errors_a2 = abs(guessesA - (-0.554924379850069));
loglog((errors_a2(1:end-1)), (errors_a2(2:end)), 'bo');
s_a2 = polyfit(log(errors_a2(1:end-1)), log((errors_a2(1:end-1))),1);
hold on
x_a2 = (errors_a2(1:end-1));
loglog( x_a2 ,  abs(s_a2(2)*10^19).*x_a2.^s_a2(1) , 'r-');
title('Convergnece Rate f(a) @ starting point = 0.15');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off

% %%% Lets try plotting each guess to see what is happening as the algorithm
% %%% progresses. Here I will loop through wach guess and plot it on the
% %%% function and we can see if and where it is getting stuck.
% figure
% %%% Original function
% plot(x_values_a, f_a(x_values_a), 'k', 'LineWidth',1.5); grid on; hold on;
% for ii=1:length(guessesA)
% xval = guessesA(ii);
% figure(5)
% plot(xval, g_a(xval), 'o', 'MarkerSize', 8); hold on;
% title(['Plotting guess ', num2str(ii)]);
% drawnow; pause(0.01);
% end
% legend('g_a(x) = -x^2 + x + 2', 'Root of g_a');
%% Plot for g_b
% figure
% plot(x_values_b, y_values_b, 'LineWidth', 1.5);
% xlabel('x');
% ylabel('g_b(x)');
% title('Function g_b(x) = exp(x) - 2 - x');

% Apply fixpt to find the root for g_b
[x_b, err,guessesB] = fixpt(f_b, -2, 1e-9, 10000); % Change 0.5 to the desired initial guess
% figure
% %%% Original function
% plot(x_values_b, f_b(x_values_b), 'r', 'LineWidth',1.5); grid on; hold on;
% 
% hold on;
% %%% What if we plot ALL of the guesses?!?
% for ii=1:length(guesses)
% figure
% plot(x_b, g_b(guesses(ii)), 'ro', 'MarkerSize', 8); hold on;
% end
% text(x_b, f_b(x_b), sprintf('Root: %.4f', x_b), 'VerticalAlignment', 'bottom');
% line([x_b, x_b], [g_b(x_b), min(y_values_b)], 'LineStyle', '--', 'Color', 'black');
% grid on;
% legend('g_b(x) = exp(x) - 2 - x', 'Root of g_b');
figure
errors_b2 = abs(guessesB - (-0.768039046747227));
loglog((errors_b2(1:end-1)), (errors_b2(2:end)), 'bo');
s_b2 = polyfit(log(errors_b2(1:end-1)), log((errors_b2(1:end-1))),1);
hold on
x_b2 = (errors_b2(1:end-1));
loglog( x_b2 ,  abs(s_b2(2)*10^15).*x_b2.^s_b2(1) , 'r-');
title('Convergnece Rate f(b) @ starting point = 0.15');
ylabel('Error(N+1)');
xlabel('Error(N)');
hold off

function y = g_a(x)
y = -x.^2+x+2;
end

function y = g_b(x)
y = exp(x) - 2 - x;
end