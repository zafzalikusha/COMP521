%% Greyscale PCA
%
% Jared Brzenski

close all;
clear all;
clc 
% Setup principal components to use. Pick 4.
principal_components = [1 10 50 100 200 400]
pics = length(principal_components);

% Read in image
I = imread('FH.jpg'); 
% and convert to grayscale
I = rgb2gray(I);

% Show it
figure(1)
imshow(I, []);
%title('A grey scaled image of Jared, His size:D 255x300 Image'); 
title('A grey scaled image-5472x2976 Image'); 

% Make data double precision
data = double(I);
%
[m n] = size(data);

% Find mean
mn   = mean(data,1);

% make data have zero mean, for covariance
data = data - repmat(mn,m,1);

fprintf('Now data has zero mean: %4.5f\n', (mean(mean(data))) );

% Find covariance of the data
covar_temp = cov(data);

% Find eigen values and vectors, returns eigenvalues as vector
[PC, V] = eig(covar_temp, 'vector');

% Flip values and vectors to go from biggest to smallest
V  = flipud(V);
PC = fliplr(PC);

% Loop over different PC sizes and plot them
pc = principal_components; 
figure(34); 
vsum = sum(V);
tiledlayout(3,2)
for pp = 1:pics
    % Extract the principal components we asked for
    output = PC(:,1:(pc(pp)))' * data';
    [xx yy] = size(output); ts = (xx+1)*yy;
    % Reconstruct full image from those principal components ( round for
    % greyscale values, need to be whole numbers )
    reconstruct = round((PC(:,1:(pc(pp)))*output) + repmat(mn,m,1)');

    % Show reconstructed image with n principal components
    nexttile
    imshow(reconstruct', []);
    title([num2str(pc(pp)) ' components, ' num2str((ts)/(m*n)*100) '%  storage, ' num2str(sum(V(1:pc(pp))/vsum)*100) '% eigen']);
end

figure(99)
    imshow(reconstruct', []);
    title([num2str(pc(pp)) ' components, ' num2str((ts)/(m*n)*100) '%  storage, ' num2str(sum(V(1:pc(pp))/vsum)*100) '% eigen']);