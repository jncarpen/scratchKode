%% January 21, 2020.
% J. Carpenter

% get a grid of locations
[X,Y] = meshgrid(1:10,1:10);
X = reshape(X, length(X)^2, 1);
Y = reshape(Y, length(Y)^2, 1);

% vector orientation
theta = 360;
thetaVec = repmat(theta, 100, 1);

% define u and v
u = cos(thetaVec .* pi/180); 
v = sin(thetaVec .* pi/180);

figure; hold on;
% set(gca, 'visible', 'off')
xlim([0 11]); ylim([0 11]);

% plot data vectors
modelVecs = quiver(X, Y, u, v);
set(modelVecs, 'Color', [.5 .7 .3]);

% set a title
tit = strcat('theta: ', sprintf('%.f', theta));
title(tit); box off;


