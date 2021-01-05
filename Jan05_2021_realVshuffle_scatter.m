% Jan 05, 2020

%% Plot real v shuffle - J.F2B
% plot the real data (entire population) against 
% the (mean) of the shuffled population

% find the mean of the shuffled distribution for each unit
clear irate_mean icont_mean
for i = 1:length(irate_array)
    irate_mean(i) = mean(irate_array{1,i}, 'omitnan');
    icont_mean(i) = mean(icont_array{1,i}, 'omitnan');
end

% find the mean of each of these distributions
mean_dist = [mean(irate_mean, 'omitnan'), mean(I_rate_content, 'omitnan')];

% calculate confidence intervals
ci_X = [mean(irate_mean) - 2*std(irate_mean), mean(irate_mean) + 2*std(irate_mean)];
ci_Y = [mean(I_rate_data) - 2*std(I_rate_data), mean(I_rate_data) + 2*std(I_rate_data)];

figure; hold on;
s = scatter(irate_mean, I_rate_data, [25], 'k', 'filled');
s.MarkerFaceAlpha = .2;
meanPoint = scatter(mean_dist(1), mean_dist(2), [60], 'g', 'filled');
meanPoint.MarkerFaceAlpha = .9;
plot(ci_X, ci_Y, '--r', 'LineWidth', 1.15)
legend('data', 'mean','95% CI'); legend boxoff
title('information rate (bits/s)', 'FontWeight','normal');
xlabel('info rate (shuffled data)'); ylabel('info rate (real data)')
set(gca,'FontSize', 15, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
box off;
