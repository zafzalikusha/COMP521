clc;
clear;
close all;
%% Define the number of terms
N = 50;
% Initialize arrays to store partial sum and values of n
partialSum = zeros(1, N);
nValues = 1:N;

% Calculate the partial sums
for n = 1:50
    partialSum(n) = sum((factorial(1:n).^2) ./ factorial(2 * (1:n)));
end
%% Plot
figure;
plot(nValues, partialSum, 'o-');
xlabel('(n)');
ylabel('(Partial Sum)');
title('Partial Sums of function (a)');
grid on;