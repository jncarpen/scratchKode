% January 05, 2020
% Autocorrelation:
% use this to show that the values we chose to shift the spiketrain by
% are logical

%% calculate stuff

[spktrn, spktrn_smooth] = binSpikes(P(:,1), ST); % get spiketrain
[xcf,lags,bounds,h] = crosscorr(spktrn_smooth, spktrn_smooth, 'NumLags',45000,'NumSTD',2);
[xcf2,lags2,~,~] = crosscorr(spktrn, spktrn, 'NumLags',45000,'NumSTD',2);
close all;

figure; set(gcf,'color','w'); hold on;
plot(lags*deltaT, xcf, 'k', 'LineWidth', 1);
plot(lags2*deltaT, xcf2, 'r', 'LineWidth', 1);
xlabel('lag (s)'); ylabel('xcf'); title('autocorrelation of spiketrain', 'FontWeight', 'normal')
set(gca,'FontSize',20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
legend('smooth', 'raw'); legend boxoff
box off;



%% plot 
% name the variables
X = I_content_shuff;
Xd = I_content;
plotName = 'information content (bits/spk)';
xname = 'info content (bits/spk)';

stdev = std(X);

% plot
figure;
set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'FaceColor', [.8 .8 .8]);
xline(Xd, '--k', 'LineWidth', 1.5);
xline(mean(X) + stdev*2, ':r', 'LineWidth', 1.5); % std
xline(mean(X) - stdev*2, ':r', 'LineWidth', 1.5); % std
ylabel("frequency (count)"); xlabel(xname);
set(gca,'FontSize',15, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
title(plotName, 'FontName', 'Helvetica UI', 'FontSize', 20, 'FontWeight', 'normal');
% l = legend('shuffled data', 'real data', '95% CI', 'Location', 'northeastoutside');
% legend boxoff    
box off;

