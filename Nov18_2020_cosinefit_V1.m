%% Run the model on V1 dataset (Wejian's data)
% November 18, 2020

%% get data
load D:\Collab\WeijianV1\WZ_V1\FullNeuronBehaviorDataSet.mat

% choose a session/unit number
sessNum = 2;
unitNum = 6;

% grab unit and corresponding [downsampled] data
unitdata = spikes_and_timestamps(FullNeuronBehaviorDataSet, sessNum, unitNum);

% store data in necessary variables
P = unitdata.position;
ST = unitdata.spiketimes;
HD = FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,sessNum}(:,9);

% visualize the cell
% figure
% pathPlot_hd(P, ST, hd)

% what is the object location (in cm)
ref_point = [32, 26]; % for session 2 (first dataset)
% ref_point = [42, 17.5]; % for session 3 (first dataset)


%% fit model
[model] = modelMe(P, ST, HD)

% calculate distance between data + prediction (in cm)
d = sqrt((ref_point(1,1)-model.bestParams.xref)^2 + (ref_point(1,2)-model.bestParams.yref)^2);

%% plot 
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');
figtit = strcat('S ', sprintf('%.f', sessNum), ', U ', sprintf('%.f', unitNum), ', Error =', sprintf('%.2f', d), ' cm');
sgtitle(figtit, 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')

subplot(1,3,1)
plot_modelDynamics(P, ST, model, ref_point)

subplot(1,3,2)
plot_iterations(model)

subplot(1,3,3)
plot_vectorMod(model)
c(2) = colorbar;
cbfreeze(jet, c(2), 'on');

fileBody = strcat('S', sprintf('%.f', sessNum), '_U', sprintf('%.f', unitNum));
filename = strcat('D:\egoAnalysis\Nov18_vectorMod_V1\', fileBody, '.png');
saveas(fig, filename);

% close all



