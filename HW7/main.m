close all;
clear;
clc;

% Set boundary conditions
c1 = 0;
c2 = 0;

% Set time step, spatial domain, and time domain parameters
dt = 0.00007812;
a = 0;
b = 1;
t0 = 0;
tf = 0.1;

% Set different ∆x values to test
dx_values = [0.2, 0.05, 0.025, 0.0125];
errors = zeros(1, length(dx_values));

% Boundary conditions
c1 = 0;
c2 = 0;

% Time-stepping loop
for dx_index = 1:length(dx_values)
    % Set the current ∆x
    dx = dx_values(dx_index);

    % Create spatial and temporal grids
    x = a:dx:b;
    t = t0:dt:tf;
    n = length(x);
    m = length(t);

    % Calculate r
    r = dt / (dx^2);

    % Initialize solution matrix
    U = zeros(n, m);

    % Initial condition
    U(:, 1) = sin(pi * x) + sin(3 * pi * x);

    % Time-stepping loop
    for j = 1:m-1
        for ii = 2:n-1
            % Explicit finite difference scheme
            U(ii, j+1) = (1 - 2*r) * U(ii, j) + r * (U(ii+1, j) + U(ii-1, j));
        end

        % Apply boundary conditions
        U(1, j+1) = 0;
        U(end, j+1) = 0;
    end

    % Calculate the exact solution at the final time step
    exact_solution = sin(pi * x) .* exp(-pi^2 * t(end)) + sin(3 * pi * x) .* exp(-9 * pi^2 * t(end));

    % Calculate the error
    errors(dx_index) = max(abs(U(:, end)' - exact_solution));

    % Create meshgrid matrices for t and x
    [T, X] = meshgrid(t, x);

    % Plot the solution using surf with logarithmic color scale
    
    figure;
    subplot(2, 1, 1);
    surf(T, X, U, 'EdgeColor', 'none');
    xlabel('Time (t)');
    ylabel('Space (x)');
    zlabel('Temperature (u)');
    title(['Numerical Solution using Explicit Finite Differences (\Deltax = ' num2str(dx) ')']);
    
    % Plot the exact solution using surf
    subplot(2, 1, 2);
    surf(T, X, sin(pi * X) .* exp(-pi^2 * t(end)) + sin(3 * pi * X) .* exp(-9 * pi^2 * t(end)), 'EdgeColor', 'none');
    xlabel('Time (t)');
    ylabel('Space (x)');
    zlabel('Exact Temperature (u)');
    title(['Exact Solution (\Deltax = ' num2str(dx) ')']);
end

% Draw the log-log plot for the error metric
figure;
loglog(dx_values, errors, 'o-', 'LineWidth', 2);
hold on;
text(dx_values, errors, arrayfun(@(e) sprintf('dx = %.4f', e), dx_values, 'UniformOutput', false), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('∆x');
ylabel('Error');
title('Error vs. ∆x');
grid on;

% Perform linear fit to find slope (order of convergence)
p = polyfit(log(dx_values), log(errors), 1);
slope = p(1);
fprintf('Slope (order of convergence): %.4f\n', slope);
