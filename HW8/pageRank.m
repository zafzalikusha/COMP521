%% Adjacency Solver

close all; clear all; clc;

A = [0, 1/2, 0, 1/2;
     1/2, 0, 1/2, 0;
     1/2, 0, 0, 0;
     0, 1/2, 1/2, 0];

d = 0.85;
[m n] = size(A);

M = d .* A + ((1-d)/n)*ones(size(A));

x = ones(n, 1);
diff = x;
count = 1; 
maxiter = 10000;

while norm(diff) > 1e-6 && count < maxiter
    xold = x;
    xnew = M * x;
    diff = xnew - xold;
    x = xnew;
    count = count + 1;
end
x_norm = x / norm(x);

disp(x_norm)

% Check with eigen values / vectors
[V L] = eig(M);

% Largest eigen value
L

% matches eigen vector (+/- sign)
V
    