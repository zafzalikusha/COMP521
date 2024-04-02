function [x,err,xx] = bisct(f,a,b,TolX,MaxIter)
%bisct.m to solve f(x) = 0 by using the bisection method.
%input : f = ftn to be given as a string ’f’ if defined in an M-file
% a/b = initial left/right point of the solution interval
% TolX = upperbound of error |x(k) - xo|
% MaxIter = maximum # of iterations
%output: x = point which the algorithm has reached
% err = (b - a)/2(half the last interval width)
% xx = history of x

TolFun=eps; fa = feval(f,a); fb = feval(f,b);
if fa*fb > 0, error('We must have f(a)f(b)<0!'); 
end

for k = 1: MaxIter
    xx(k) = (a + b)/2;
    fx = feval(f,xx(k)); err = (b-a)/2;
    if abs(fx) < TolFun | abs(err)<TolX, break;
    elseif fx*fa > 0, a = xx(k); fa = fx;
    else b = xx(k);
    end
end

x = xx(k);
if k == MaxIter, fprintf('The best in %d iterations\n',MaxIter),
end