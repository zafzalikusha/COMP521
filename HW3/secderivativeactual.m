% File: secderivativeactual.m
function Dactual = secderivativeactual(xgrid)
% Function that evaluates the exact second derivative of a function f(x)
% on a set of grid points
%
% Inputs:
% xgrid : x points at which the second derivative is evaluated
%
% Outputs:
% Dactual : second derivative values at the grid points

Dactual = 100 * cos(10 * xgrid) + 2;
end
