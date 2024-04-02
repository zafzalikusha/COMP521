function [x,err,xx] = fixpt(g,x0,TolX,MaxIter)

% solve x = g(x) starting from x0 by fixed-point iteration.
%input : g,x0 = the function and the initial guess
% TolX = upperbound of incremental difference |x(n + 1) - x(n)|
% MaxIter = maximum # of iterations
% output: x = point which the algorithm has reached
% err = last value |x(k) - x(k - 1)| achieved
% xx = history of x

if nargin < 4, MaxIter = 100; end
if nargin < 3, TolX = 1e-9; end
xx(1) = x0;

for k = 2:MaxIter
    xx(k) = feval(g,xx(k - 1)); %Eq.(4.1.3)
    err = abs(xx(k) - xx(k - 1)); if err < TolX, break; end
end
x = xx(k);

if k == MaxIter
    fprintf('Do not rely on me, though best in %d iterations\n',MaxIter)
end