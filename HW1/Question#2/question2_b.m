clc;
clear;
close all;
%% The symbolic variable and the function
syms x;
f = log(3*x + 1);
% Create an array of x values, x is between 0 and 0.5, with delta x=0.001
x_values = 0:0.001:0.5;
%% Initialize arrays to store error values
error_3_terms = zeros(size(x_values));
error_6_terms = zeros(size(x_values));
% Compute the error for each x value
for i = 1:length(x_values)
    % Compute the exact value
    exact_value = double(subs(f, x, x_values(i)));

    % Compute the Taylor series approximations with 2 and 5 terms
    taylor_3 = taylor(f, x, 'Order', 3);
    taylor_6 = taylor(f, x, 'Order', 6);
    
    % Evaluate the approximations at the current x value
    approx_3 = double(subs(taylor_3, x, x_values(i)));
    approx_6 = double(subs(taylor_6, x, x_values(i)));
    
    % Compute the error and store it in the arrays
    error_3_terms(i) = abs(exact_value - approx_3);
    error_6_terms(i) = abs(exact_value - approx_6);
end

% Compute the L2-norm of the error vectors
norm_2_terms = sqrt(sum(error_3_terms.^2));
norm_5_terms = sqrt(sum(error_6_terms.^2));
%% Display the results
fprintf('L2-Norm of Error (2 Terms): %.6f\n', norm_2_terms);
fprintf('L2-Norm of Error (5 Terms): %.6f\n', norm_5_terms);