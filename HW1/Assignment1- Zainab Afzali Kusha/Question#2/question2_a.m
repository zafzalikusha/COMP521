clc;
clear;
close all;
%% The symbolic variable and the function
syms x;
f = log(3*x + 1);
% Taylor series with 2, and 5 terms
taylor_2 = taylor(f, x, 'Order', 2); 
taylor_3 = taylor(f, x, 'Order', 3); 
taylor_5 = taylor(f, x, 'Order', 5); 
taylor_6 = taylor(f, x, 'Order', 6); 
%% Convert the symbolic expressions to functions
taylor_3_func = matlabFunction(taylor_3);
taylor_6_func = matlabFunction(taylor_6);
% Create an array of x values, x is between 0 and 0.5, with delta x=0.001
x_values = 0:0.001:0.5;
% Evaluate the Taylor series functions at the x values
y_taylor_3 = taylor_3_func(x_values);
y_taylor_6 = taylor_6_func(x_values);
% Calculate the exact function values
y_exact = log(3*x_values + 1);
%% Plot the results
figure;
plot(x_values, y_exact, 'b', 'LineWidth', 2, 'DisplayName', 'Exact Function');
hold on;
plot(x_values, y_taylor_3, '--g', 'LineWidth', 2, 'DisplayName', 'Taylor Series (2 terms)');
plot(x_values, y_taylor_6, '--r', 'LineWidth', 2, 'DisplayName', 'Taylor Series (5 terms)');
xlabel('(x)');
ylabel('(y)');
title('Exact Function and Taylor Series Approximations');
legend('Location', 'NorthWest');
grid on;
hold off;