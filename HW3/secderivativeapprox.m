function [xgrid, Dapprox, aproxlim] = secderivativeapprox(x, h, Fx)
% This function implements the finite difference approximation of the
% second derivative of a function f(x)
%
% Inputs:
% x   : the X interval
% h   : the step size
% Fx  : handle to the function f(x)
%
% Outputs:
% xgrid     : The grid points in X
% Dapprox   : The finite difference approximation at the grid points
% aproxlim  : The indices of the interval endpoints where the approximation
%             is calculated

% Create the vector with the grid points
xgrid = x(1):h:x(2);

% Set the approximation limits. These limits depend on the finite
% difference approximation you will use.
% Note: In this example, we cannot include the first and last endpoints
% into the approximation
% firstendpoint variable represents the index of the first grid point where the calculation of the second derivative begins.
% In the code, it is set to 3, which means that the calculation starts from the third grid point in the xgrid array. 
% Setting it to 3 ensures that the calculation avoids the first two grid points because the 
% finite difference method requires neighboring points for the approximation, 
% and those points would not have enough neighbors.
% lastendpoint variable represents the index of the last grid point where the calculation of the second derivative ends.
% It is set to size(xgrid, 2) - 2. The size(xgrid, 2) expression gives the total number of grid points in the xgrid array.
% Subtracting 2 ensures that the calculation stops at the second-to-last grid point.
% Like the first point, this is done to ensure that there are enough neighboring points
% for the finite difference calculation and to avoid reaching beyond the last point.

firstendpoint = 3;
lastendpoint = length(xgrid)-2;

% Initialize the array approxlim
aproxlim = [3, length(xgrid)-2];

% Initialize the vector that will have the approximation 
% and be the exact size of xgrid
Dapprox = zeros(size(xgrid));

% Calculate the function values
% this loop iterates through the elements of xgrid, evaluates the function Fx at each xgrid value, 
% and stores the results in the Dapprox array. 
    for i = firstendpoint:lastendpoint
     Dapprox(i) = Fx(xgrid(i));
    end
% Calculate the finite difference approximation
    dx = 1 / (12 * h^2);

    for i = 3:length(xgrid)-2
    Dapprox(i) = (-Fx(xgrid(i + 2)) + 16 * Fx(xgrid(i + 1)) - 30 * Fx(xgrid(i)) + 16 * Fx(xgrid(i - 1)) - Fx(xgrid(i - 2))) * dx;
    end

end
