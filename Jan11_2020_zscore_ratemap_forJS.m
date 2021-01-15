% Jan 11, 2020
% to show JS what the z-scored ratemap deal is all about

% grab some position data
whichSession = randi(length(JZ.neurons)); % randomly choose a session
P_raw = units2cm(JZ.neurons(whichSession).members(1).P);
P_smooth = smooth_pos(P_raw, 2);

% set parameters
clear simparam
simparam.position = P_smooth;
simparam.ctr_mass = [rand(1)*150, rand(1)*150];
simparam.ref_point = [rand(1)*150, rand(1)*150];
simparam.noise = 0; % noiseLevel(simcell);
simparam.width = rand(1)*5 + 5;
simparam.theta = rand(1)*360;

% simulate cell
[sim] = simulate_place(simparam);

% set spiketimes
ST = sim.spiketimes;

% sampling frequency (50 samples/sec)
Fs = mode(diff(P_raw(:,1))); 

% get 'real' head direction values (deg)
HD = get_hd(P_smooth);

% make a ratemap
figure; set(gcf,'color','w');
map = analyses.map(P_smooth, ST, 'smooth', 2, 'binWidth', 150/75); % calculate tuning curve
newMap = (map.z-M)./V;
newMap = (newMap-mean(squeeze(newMap), 'omitnan'))./std(newMap, [0], 'all', 'omitnan');
peakRate = nanmax(nanmax(newMap));
rate_map_title = strcat('max zscore: ', sprintf('%.2f',peakRate));
plot.colorMap(newMap);
pbaspect([1 1 1])
c = colorbar; c.FontName = 'Helvetica'; c.FontSize = 15;
colormap(gca,'jet')
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Helvetica', 'FontSize', 15, 'FontWeight', 'normal');
box off

% find mean
M = mean(map.z, 'all', 'omitnan');

% find variance 
V = var(map.z,[],'all', 'omitnan');


