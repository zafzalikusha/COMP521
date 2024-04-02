function x = solveLinearEquations(b, A, indxi)
    % Solves a system of linear equations A * x = b using LU decomposition
    % with partial pivoting.
    % Inputs:
    %   b: Right-hand side vector (N x 1)
    %   A: LU decomposed matrix (N x N)
    %   indxi: Permutation vector (N x 1)
    % Output:
    %   x: Solution vector (N x 1)

    N = length(b);
    x = zeros(N, 1);
    
    % Forward substitution
    % Initialize ii
    ii = 0; 
    for i = 1:N
        ip = indxi(i);
        sum = b(ip);
        b(ip) = b(i);
        
        if ii > 0
            for j = ii:i - 1
                sum = sum - (A(i, j) * b(j));
            end
        elseif sum ~= 0
            ii = i;
        end
        
        b(i) = sum;
    end
    
    % Backward substitution
    for i = N:-1:1
        sum = b(i);
        for j = (i + 1):N
            sum = sum - A(i, j) * x(j);
        end
        x(i) = sum / A(i, i);
    end
end
