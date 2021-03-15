% load('D:\Data\External Data\Simulation-EX\ego0.mat')
% load('D:\Data\External Data\Simulation-EX\ego90.mat')
% load('D:\Data\External Data\Simulation-EX\ego180.mat')
% load('D:\Data\External Data\Simulation-EX\ego270.mat')
% load('D:\Data\External Data\Simulation-EX\hd0.mat')
% load('D:\Data\External Data\Simulation-EX\hd90.mat')
% load('D:\Data\External Data\Simulation-EX\hd180.mat')
% load('D:\Data\External Data\Simulation-EX\hd270.mat')
% load('D:\Data\External Data\Simulation-EX\place.mat')
% load('D:\Data\External Data\Simulation-EX\placeego90.mat')
% load('D:\Data\External Data\Simulation-EX\placehd270.mat')
% load('D:\Data\External Data\Simulation-EX\egodist270.mat')

% RP = [75 75];
pathPlot_hd(P, ST, get_hd(P));
out = modelMe(P,ST,get_hd(P));
plotMe(out);

outRM = plot_egoRM(P, ST, RP);
ego_RM = egoRateMap(P, ST, RP);
d = sqrt((0 - (10)).^2 + (0 - 10).^2);


plot_MVLmap(P,ST,HD)

tc_shuffle(P, ST, RP)

figure
map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 150/75); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z); shading interp
pbaspect([1 1 1])
% colormap(gca,'jet')
set(gca, 'visible', 'off', 'box', 'off');
% c2 = colorbar; c2.FontSize = 25;
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
box off

[strain] = binSpikes(P(:,1), ST);
hdTuning(HD, P, strain);






