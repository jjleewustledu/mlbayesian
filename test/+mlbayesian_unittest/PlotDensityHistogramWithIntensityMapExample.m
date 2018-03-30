%% Plot Density Histogram with Intensity Map
% Load the sample data. 

% Copyright 2015 The MathWorks, Inc.

load seamount 
%% 
% Correct grid for negative y-values and draw histogram in 2D.
hold on
dat = [-y,x];  
hist3(dat)    
%%
% Extract histogram data.
n = hist3(dat); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
%%
% Generate grid for 2-D projected view of intensities.
xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
%%
% Make a pseudocolor plot.
h = pcolor(xb,yb,n1);
%%
% Set the z-level and colormap of the displayed grid, and display the
% default 3-D perspective view.
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot) % heat map 
title('Seamount:Data Point Density Histogram and Intensity Map');
grid on 
view(3);