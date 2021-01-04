%% Simulate 100 place cells and run the model on them

P = data(4).pos;
step = 0.5;
x = min(P(:,2))+10:step:max(P(:,2))-10;
y = min(P(:,3))+10:step:max(P(:,3))-10;

clear varEx tuningStrength pc100 egoCell_100

for i = 1:100
clear param sim model
% set parameters for simulated cell
param.position = data(4).pos;
x_i = randsample(x,1)';
y_i = randsample(y,1)';
param.ctr_mass = [x_i y_i]; % center of arena
param.ref_point = [randsample(x,1) randsample(y,1)];
param.noise = 0; % no noise 
param.width = 1;
param.theta = 270; %rand(1).*360;
param.radius = 40; %rand(1).*80;
%simulate the cell
[sim] = simulate_place(param);
% run the model
[model] = modelMe(sim.position, sim.spiketimes);
ref_point = [model.bestParams.xref, model.bestParams.yref]
model.bestParams.thetaP
figure
plot_vectorMod_model(model)

%% save stuff
varEx(i) = model.varExplained.model;
tuningStrength(i) = model.bestParams.g;
egoCell_100(i).st = sim.spiketimes;
egoCell_100(i).model = model;
egoCell_100(i).param = param;

% visualize
subplot(1,3,1)
plot_modelDynamics(sim.position, sim.spiketimes, model, ref_point);
subplot(1,3,2)
plot_vectorMod(model)
subplot(1,3,3)
plot_iterations(model)

% get spiketrain
[spikeTrain, ~] = binSpikes(sim.position(:,1), sim.spiketimes);
[binCtrs, tcVals] = goalDirSar(sim.position, ref_point, get_hd(sim.position), sim.spiketimes, 40);

fig = figure('units','normalized','outerposition',[0 0 1 1]); % make fullscreen fig
set(gcf,'color','w');
fileBody = strcat('allo_', sprintf('%.f', i));
subplot(2,3,1)
plot_modelDynamics(sim.position, sim.spiketimes, model, ref_point);
hold on; plot(73, 112, 'o')
subplot(2,3,2)
plot(binCtrs, tcVals, 'LineWidth', 1.1, 'Color', 'k'); box off;
xlabel('egocentric bearing (deg)'); ylabel('firing rate (Hz)')
subplot(2,3,3)
plot_iterations(model)
subplot(2,3,4)
plot_vectorMod_data(model)
subplot(2,3,5)
plot_vectorMod(model)
subplot(2,3,6)
plot_vectorMod_model(model)

% save figures
filename = strcat('D:\egoAnalysis\test\', fileBody, '.png');
saveas(fig, filename);

% click through all figures
close all 
end

model = JZM.neurons(105).members(1).model;


% make
figure
hold on;
myBins = -3:.25:3;
% myBins = -.5:.05:.5
histogram(TS.place, myBins, 'FaceColor', 'b');   
ax = gca; alpha(ax,.1);
histogram(TS.ego, myBins, 'FaceColor', 'r');
histogram(TS.egoDist, myBins, 'FaceColor', 'g');
histogram(TS.hd, myBins, 'FaceColor', 'y');
histogram(TS.placeMod, myBins, 'FaceColor', 'm');
histogram(TS.placeEgo, myBins, 'FaceColor', 'c');
title("Model-Predicted Tuning Strength", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
xlabel("tuning strength (g)", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
ylabel("frequency", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
legend('place', 'ego bearing', 'ego bearing + dist', 'head direction', 'placeMod','placeEgo')



%% VARIANCE EXPLAINED
hold on;
myBins = 0:.05:1;
histogram(VE.place, myBins, 'FaceColor', 'b');   
ax = gca; alpha(ax,.1);
histogram(VE.ego, myBins, 'FaceColor', 'r');
histogram(VE.egoDist, myBins, 'FaceColor', 'g');
histogram(VE.hd, myBins, 'FaceColor', 'y');
histogram(VE.placeMod, myBins, 'FaceColor', 'm');
histogram(VE.placeEgo, myBins, 'FaceColor', 'c');
title("Variance Explained by Model", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
xlabel("variance explained", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
ylabel("frequency", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
% legend('place', 'hd-mod place', 'ego', 'hd', 'ego+dist', 'allo', 'FontSize', 14)
legend('place', 'ego bearing', 'ego bearing + dist', 'head direction', 'placeMod','placeEgo')
xticks([0 .2 .4 .6 .8 1])

%% egocentric v. allocentric cells
egocells = [VE.egoDist, VE.ego, VE.placeEgo];
allocells = [VE.place, VE.hd, VE.allo];
myBins = -3:.25:3;
histogram(TS.place, myBins, 'FaceColor', 'k');
hold on;
ax = gca; alpha(ax,.1);
histogram(TS.ego, myBins, 'FaceColor', 'r');
legend('pure place cells', 'egocentric bearing cells')





% make plots
pathPlot_hd(sim.position, sim.spiketimes, get_hd(sim.position))
pathPlot_hd(P, ST, get_hd(P))
hold on; plot(param.ref_point(1,1), param.ref_point(1,2), 'o')
% hold on; plot(103, 82, 'o')


pathPlot_quiver(sim.position, sim.spiketimes, get_hd(sim.position))
pathPlot_quiver(P, sim.spiketimes, get_hd(sim.position))

figure
tc_shuffle(sim.position, sim.spiketimes, ref_point)

figure
tc_shuffle_poisson(sim.position, sim.spiketimes, ref_point, 100);

figure 
ref = [0 0]
egoRateMap(sim.position, sim.spiketimes, ref_point)


figure
map = analyses.map(sim.position, sim.spiketimes, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
colormap(gca,'jet')
c2 = colorbar; c2.FontSize = 25;
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
box off



%% distance formula

modelx = 427680;
modely = 207780;
simx = 64;
simy = 131;

d = sqrt((simx-modelx)^2 + (simy-modely)^2)

