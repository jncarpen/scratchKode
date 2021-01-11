% Jan 05, 2020

%% Plot real v shuffle - J.F2B
% plot the real data (entire population) against 
% the (mean) of the shuffled population

X_data = MS_HD;
X_shuf = MS_HD_shuff;
plotTitle = 'RH Modulation (Real v. Shuffle)';
xname = 'rh mod (shuffle)';
yname = 'rh mod (data)';

% find the mean of each of these distributions
mean_dist = [mean(X_shuf, 'omitnan'), mean(X_data, 'omitnan')];

% calculate confidence intervals
ci_shuf = [mean(X_shuf) - 2*std(X_shuf), mean(X_shuf) + 2*std(X_shuf)];
ci_data = [mean(X_data) - 2*std(X_data), mean(X_data) + 2*std(X_data)];

figure; hold on;
s = scatter(X_shuf, X_data, [25], 'k', 'filled');
s.MarkerFaceAlpha = .2;
meanPoint = scatter(mean_dist(1), mean_dist(2), [60], [0 .8 .2], 'filled');
meanPoint.MarkerFaceAlpha = .6;
plot(ci_shuf, ci_data, '--r', 'LineWidth', 1.15) % CI line
plot(linspace(0, 6), linspace(0, 6), ':b', 'LineWidth', 1.15); % equality line
legend('data', 'mean','95% CI'); legend boxoff
title(plotTitle, 'FontWeight','normal');
xlabel(xname); ylabel(yname); xlim([0 1]); ylim([0 1]);
set(gca,'FontSize', 20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
box off;
