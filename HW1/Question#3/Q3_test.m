%test 1 and test 2
b1 = [2*pi; 5*pi; -8*pi]
b2 = [3 ; 0; -1]

% Example matrix A (already LU decomposed)
A1 = [1 2 -1; 6 -5 4; -9 8 -7]
A2 = [pi 3*pi 2*pi; 0 -1 2/3; -pi -3*pi 2*pi]

% Permutation vector indxi ???
% Perform LU decomposition with partial pivoting
[L, U, P1] = lu(A1);
[L, U, P2] = lu(A2);

% Get the permutation vector indxi
indxi_1 = sum(P1, 2);
indxi_2 = sum(P2, 2);

% Solve the linear equations using forward-backward substitution
x1 = solveLinearEquations(b1, A1, indxi_1);
x2 = solveLinearEquations(b2, A2, indxi_2);

disp('MATLAB Solution Vector x1:');
disp(A1\b1); disp('MATLAB Solution Vector x2:');
disp(A1\b2); 