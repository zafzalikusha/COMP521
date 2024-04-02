function [indxi, d, A] = luDecomposition(A)
    % LU decomposition with rowwise partial pivoting
    % Inputs:
    %   A: Input matrix (N x N)
    % Outputs:
    %   A: LU decomposed matrix (N x N)
    %   indxi: Permutation vector (N x 1)
    %   d: Determinant (+1 or -1)
    % Get the size of the matrix a
    N = size(A, 1);

    % Initialize variables
    indxi = 1:N;
    d = 1;
    TINY = 1.0e-20;

    % Initialize vv as a vector of ones
    vv = ones(N, 1); %???

    for i = 1:N
        big = 0.0;

        for j = 1:N
            temp = abs(A(i, j));
            if temp > big
                big = temp;
            end
        end

        if big == 0.0
            error('Singular matrix');
        end

        vv(i) = 1.0 / big;
    end

    for j = 1:N
        for i = 1:j-1
            sum = A(i, j);
            for k = 1:i-1
                sum = sum - (A(i, k) * A(k, j));
            end
            A(i, j) = sum;
        end

        big = 0.0;

        for i = j:N
            sum = A(i, j);
            for k = 1:j-1
                sum = sum - (A(i, k) * A(k, j));
            end
            A(i, j) = sum;

            dum = vv(i) * abs(A(i, j));

            if dum >= big
                big = dum;
                imax = i;
            end
        end

        if j ~= imax
            for k = 1:N
                dum = A(imax, k);
                A(imax, k) = A(j, k);
                A(j, k) = dum;
            end

            d = -d;
            vv(imax) = vv(j);
        end
            indxi(j) = imax; %??
        
        if abs(A(j, j)) == 0
            A(j, j) = TINY;
        end

        if j ~= N
            dum = 1.0 / A(j, j);
            for i = j+1:N
                A(i, j) = A(i, j) * dum;
            end
        end
    end
end
