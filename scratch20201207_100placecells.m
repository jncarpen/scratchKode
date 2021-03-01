%% Simulate 100 place cells and run the model on them

for i = 1:100
% root structure for ratemap
P = dataof(i).P;
root.A = 10;
root.ctr = [75 75];
root.sigma = [14 18];
root.size = 150;
root.bins = 20;
root.P = P;

%simulate the cell
[map] = simulate_ratemap(root);
[sim] = simulate_place(map,P(:,1));

% run the model
[out] = modelMe(P, sim.ST, get_hd(P));

% save
place100new(i).root = root;
place100new(i).ST = sim.ST;
place100new(i).out = out;
end


%% Simulate 100 egocentric cells and run the model on them
for i = 1:100
% parameters for simulation
P = dataof(i).P;
param.theta = randi(360);
param.P = P;
param.Z = get_hd(P);
param.kappa = 5;
param.rp = [75 75];
param.A = 8;

%simulate the cell
[sim, ~] = simulate_ego(param);

% run the model
[out] = modelMe(P, sim.ST, get_hd(P));

% save
ego100new(i).root = param;
ego100new(i).ST = sim.ST;
ego100new(i).out = out;
end

%% save stuff
varEx(i) = out.varExplained.model;
tuningStrength(i) = out.bestParams.g;
egoCell_100(i).st = sim.spiketimes;
egoCell_100(i).model = out;
egoCell_100(i).param = param;

% visualize
subplot(1,3,1)
plot_modelDynamics(sim.position, sim.spiketimes, out, ref_point);
subplot(1,3,2)
plot_vectorMod(out)
subplot(1,3,3)
plot_iterations(out)

% get spiketrain
[spikeTrain, ~] = binSpikes(sim.position(:,1), sim.spiketimes);
[binCtrs, tcVals] = goalDirSar(sim.position, ref_point, get_hd(sim.position), sim.spiketimes, 40);

fig = figure('units','normalized','outerposition',[0 0 1 1]); % make fullscreen fig
set(gcf,'color','w');
fileBody = strcat('allo_', sprintf('%.f', i));
subplot(2,3,1)
plot_modelDynamics(sim.position, sim.spiketimes, out, ref_point);
hold on; plot(73, 112, 'o')
subplot(2,3,2)
plot(binCtrs, tcVals, 'LineWidth', 1.1, 'Color', 'k'); box off;
xlabel('egocentric bearing (deg)'); ylabel('firing rate (Hz)')
subplot(2,3,3)
plot_iterations(out)
subplot(2,3,4)
plot_vectorMod_data(out)
subplot(2,3,5)
plot_vectorMod(out)
subplot(2,3,6)
plot_vectorMod_model(out)

% save figures
filename = strcat('D:\egoAnalysis\test\', fileBody, '.png');
saveas(fig, filename);

% click through all figures
close all 
end

out = JZM.neurons(105).members(1).model;


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

