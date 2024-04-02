clear;
close all;
clc;
format long;

N = 81;
h_values = [0.01, 0.001, 0.0001];
x1 = 0.0;
x2 = 1.0;

num_iterations = 1;  % Number of iterations to run. For calculating the mean time, I used 10 iterations.

% Create cell arrays to store errors and times for each method
errors_first_derivative = cell(length(h_values), 5);
errors_second_derivative = cell(length(h_values), 5);
times_first_derivative = cell(length(h_values), 5, num_iterations);
times_second_derivative = cell(length(h_values), 5, num_iterations);

% Loop over different step sizes
for idx_h = 1:length(h_values)
    h = h_values(idx_h);
    
    % Loop over different accuracy orders (O)
    for O = 1:5
        x = x1:h:x2;
        f = exp(-5.*x).*sin(75.*x);
        f = f(:);
        x = x(:);
        % Generate random noise
        noise = 0.0001 * rand(size(x));

        % Add noise to the original function
        f_noisy = f + noise;

        % Calculate exact derivatives
        f_prime_exact_values = -exp(-5.*x).*(5.*sin(75.*x)-75*cos(75.*x));
        f_double_prime_exact_values = -exp(-5.*x).*(5600*sin(75.*x)+750*cos(75.*x));

        % Run the computations multiple times
        for iter = 1:num_iterations
            % Estimate first derivative and measure time
            tic;
            f_prime_estimate = f_derivative(1, O, f_noisy, h);
            times_first_derivative{idx_h, O, iter} = toc;

            % Estimate second derivative and measure time
            tic;
            f_double_prime_estimate = f_derivative(2, O, f_noisy, h);
            times_second_derivative{idx_h, O, iter} = toc;
        end

        % Calculate errors
        error_first_derivative = abs(f_prime_exact_values- f_prime_estimate);
        error_second_derivative = abs(f_double_prime_exact_values - f_double_prime_estimate);

        % Store errors in cell arrays
        errors_first_derivative{idx_h, O} = error_first_derivative;
        errors_second_derivative{idx_h, O} = error_second_derivative;
    end
end

% Plot the errors
figure;

% Plot for the first derivative
subplot(2, 1, 1);
for O = 1:5
    mean_errors = cellfun(@(x) mean(x), errors_first_derivative(:, O));
    loglog(h_values, mean_errors, '-o', 'DisplayName', ['Order ', num2str(O)]);
    hold on;
end
title('Mean Errors for First Derivative');
xlabel('Step Size (h)');
ylabel('Mean Error');
legend('Location', 'Best');
grid on;

% Plot for the second derivative
subplot(2, 1, 2);
for O = 1:5
    mean_errors = cellfun(@(x) mean(x), errors_second_derivative(:, O));
    loglog(h_values, mean_errors, '-o', 'DisplayName', ['Order ', num2str(O)]);
    hold on;
end
title('Mean Errors for Second Derivative');
xlabel('Step Size (h)');
ylabel('Mean Error');
legend('Location', 'Best');
grid on;

mean_times_first_derivative = zeros(length(h_values), 5);
mean_times_second_derivative = zeros(length(h_values), 5);

for idx_h = 1:length(h_values)
    for O = 1:5
       mean_times_first_derivative(idx_h, O) = mean(cell2mat(times_first_derivative(idx_h, O, :)));
       mean_times_second_derivative(idx_h, O) = mean(cell2mat(times_second_derivative(idx_h, O, :)));
    end
end

% Plot mean times for first derivative
figure;
subplot(2, 1, 1);
for O = 1:5
    loglog(h_values, mean_times_first_derivative(:, O), '-o', 'DisplayName', ['Order ', num2str(O)]);
    hold on;
end
title('Mean Times for First Derivative');
xlabel('Step Size (h)');
ylabel('Mean Time (seconds)');
legend('Location', 'Best');
grid on;

% Plot mean times for second derivative
subplot(2, 1, 2);
for O = 1:5
    loglog(h_values, mean_times_second_derivative(:, O), '-o', 'DisplayName', ['Order ', num2str(O)]);
    hold on;
end
title('Mean Times for Second Derivative');
xlabel('Step Size (h)');
ylabel('Mean Time (seconds)');
legend('Location', 'Best');
grid on;
%%
function Diff_fun = f_derivative(m, O, f, h)

n = m + O;
C = Coeff(m, O, h);
Bickley_scale = h^m * factorial(n-1) / factorial(m);
Bickley_coeff = C * Bickley_scale;
Bickley_coeff = int64(Bickley_coeff);

% Modifying the numerical stencil for a unified accuracy
if (mod(m, 2) == 0 && mod(O, 2) ~= 0)
    np = n + 1;
    j = np / 2;
    for i = 1:j - 1
        coef_u(i, :) = [C(i, :) 0];
    end

    coef_m(1, :) = [0 C(j - 1, :)];
    coef_m(2, :) = [C(j + 1, :) 0];

    for i = j + 2:np
        coef_l(i - j - 1, :) = [0 C(i - 1, :)];
    end

    C = [coef_u; coef_m; coef_l];

else
    np = n;
end

% Generating the derivative over the real domain of N nodes
[N b] = size(f);

if N < np
    error('The size of the given data should be greater than or equal to %d. \n\t Try using larger set of data or decrease the accuracy order.', np);
end

if mod(np, 2) == 0
    j = np / 2;
else
    j = (np + 1) / 2;
end

U = [C(1:j - 1, :) zeros(j - 1, N - np)];
L = [zeros(np - j, N - np) C(j + 1:np, :)];
k = C(j, :);
M = zeros(N - np + 1, N);

for i = 1:N - np + 1
    for j = 1:N
        if i == j
            M(i, j:j + np - 1) = k;
        end
    end
end

coef = [U; M; L];
Diff_fun = coef * f;
end

function C = Coeff(m, O, h)

n = m + O;
k = [1:n];
k = k(:);

for i = 1:n
    v = k - i;
    sigma = Sig(v, n);

    for l = 1:n
        C(i, l) = ((-1)^(l - m - 1) * factorial(m) / (factorial(l - 1) * factorial(n - l) * h^m)) * sigma(n - m, l);
    end
end
end

function S = Sig(v, n)

S = ones(n, n);
s = S;

for k = 1:n
    vv = v;
    vv(k) = v(1);

    for i = 2:n
        s(i, i) = s(i - 1, i - 1) * vv(i);

        if i > 2
            for j = i - 1:-1:2
                s(j, i) = s(j - 1, i - 1) * vv(i) + s(j, i - 1);
            end
        end
    end

S(:, k) = s(:, n);
end
end