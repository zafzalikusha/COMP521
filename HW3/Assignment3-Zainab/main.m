%% MAIN Program for HW 3
% COMP 521
%
% 1) Calculate the finite difference approximation for the second
%    derivative of a function f(x) on the interval x \in [ 0.4 , 1]
% 2) Determine the order of accuracy of the finite difference
%    approximation using:
%
%   2.1) Absolute error of the midpoint
%   2.2) The root mean square error of the approximation
%   2.3) The infinity norm
%
%   You have to plot the error metrics versus 4 different grid sizes for
%   each case 2.1, 2.2, 2.3. Use loglog plot. Fit the points in the loglog
%   plot to a straight line using polyfit.
%
close all; clear; clc; 

% Use the following grid sizes
h = [ 0.1 ; 0.05 ; 0.025 ; 0.0125 ];

% Calculate the number of grid sizes
m = size( h , 1 );

% Specify the Interval
x = [ 0.4 ; 1];

% Initialize an array with the error metrics
errorh = zeros( m , 3 );

% Calculate teh approzimation error for different grid sizes

for i = 1:m
    
    % Apply finite difference approximation
    [ xgrid, Dapprox, aproxlim ] = secderivativeapprox( x, h(i), @Fx );
    
    % Calculate exact solution at the ggrid points
    Dactual = secderivativeactual( xgrid );
    
    % Calculate the error vector
    % Note: Only inlcude the grid points in which tthe approximation was 
    %       computed

    % Calculate the error vector
    % Exclude the first two and last two points of xgrid
    Error = abs(Dapprox(3:end-2) - Dactual(3:end-2));
    % For the midpoint INDEX!!
    impoint = round( ( x(2) - x(1) ) / ( 2*h(i))) -1;
    
    % Show the midpoint you are using for the current grid to verify
    % that you are using the same x at each h
    fprintf('Midpoint for h=%11.10f is x=%11.10f \n',h(i),xgrid(impoint));
    
    % Note: Exclude the first two and last two points
    % Calculate the error vector
    errorh(i,1) = Error(impoint);
    errorh(i,2) = sqrt(sum(Error.^2) / numel(Error));
    errorh(i,3) = max(Error);

end
%% Plot your results
tiledlayout(3,1);

nexttile
loglog(h, errorh(:,1), '-o', 'LineWidth', 1.5, 'MarkerSize', 8); % Customize line and marker attributes
grid on; % Display grid lines
xlabel('Grid size (h)');
ylabel('|Error_{mid}|');
title('Midpoint Error');

nexttile
loglog(h, errorh(:,2), '-s', 'LineWidth', 1.5, 'MarkerSize', 8); % Customize line and marker attributes
grid on;
xlabel('Grid size (h)');
ylabel('RMSE Error');
title('Root Mean Square Error');

nexttile
loglog(h, errorh(:,3), '-d', 'LineWidth', 1.5, 'MarkerSize', 8); % Customize line and marker attributes
grid on;
xlabel('Grid size (h)');
ylabel('Infinity Norm Error');
title('Infinity Norm Error');

% Add a title to the entire figure
sgtitle('Accuracy Analysis of Finite Difference Approximation');
%% Verify with linear plot fitting
disp(' ' );
% Midpoint Error
Efit = polyfit( log(h), log(errorh(:,1)),1);
disp(['Midpoint: Fit is |E| = ' num2str(Efit(1)) '*h + (' num2str(Efit(2)) ')' ]);
% RMSE Error
Efit = polyfit( log(h), log(errorh(:,2)), 1);
disp(['RMSE: Fit is |E| = ' num2str(Efit(1)) '*h + (' num2str(Efit(2)) ')' ]);
% Infinity Error
Efit = polyfit( log(h), log(errorh(:,3)), 1);
disp(['Infinity: Fit is |E| = ' num2str(Efit(1)) '*h + (' num2str(Efit(2)) ')' ]);