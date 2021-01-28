
% load
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example1_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example2_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example3_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example4_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example5_RawData.mat')

% cell 1
figure; hold on;
plot(Cell1_Trajectory(:,2), Cell1_Trajectory(:,3), 'Color', [.8 .8 .8]);
scatter(Cell1_Spikes(:,2), Cell1_Spikes(:,3), '.');

% cell 2
figure; hold on;
plot(Cell2_Trajectory(:,2), Cell2_Trajectory(:,3), 'Color', [.8 .8 .8]);
scatter(Cell2_Spikes(:,2), Cell2_Spikes(:,3), '.');

% cell 3
figure; hold on;
plot(Cell3_Trajectory(:,2), Cell3_Trajectory(:,3), 'Color', [.8 .8 .8]);
scatter(Cell3_Spikes(:,2), Cell3_Spikes(:,3), '.');

% cell 4
figure; hold on;
plot(Cell4_Trajectory(:,2), Cell4_Trajectory(:,3), 'Color', [.8 .8 .8]);
scatter(Cell4_Spikes(:,2), Cell4_Spikes(:,3), '.');

% cell 5
figure; hold on;
plot(Cell5_Trajectory(:,2), Cell5_Trajectory(:,3), 'Color', [.8 .8 .8]);
scatter(Cell5_Spikes(:,2), Cell5_Spikes(:,3), '.');

% get movement direction
trajNow = Cell1_Trajectory;
t = trajNow(:,1);
y1 = trajNow(:,3);
x1 = trajNow(:,2);

for ts = 1:length(t)-1  
    if ts == 1
        MD(1) = NaN;
    elseif ts == length(t)-1
        MD(end-1) = NaN;
        MD(end) = NaN;
        MD(end+1) = NaN;
    else
        MD(ts) = atan2d( y1(ts+1) - y1(ts-1), x1(ts+1) - x1(ts-1))+180;
    end
end

% choose a neuron
P = Cell1_Trajectory;
ST = Cell1_Spikes(:,1);

% figure
% pathPlot_quiver(P, ST, MD)

% savefig(fig, 'D:\egoAnalysis\Jan20_PabloCells\Cell5.fig')

% run the model
total_iters = 100; 
clear monte error_for_comparison
for optim_iter = 1:total_iters 
    % choose some initial conditions randomly
    initial = choose_initial_conditions(P);
    % run the model 
    [model] = modelMe(P, ST, MD, initial);

    monte(optim_iter).model = model; % save all the runs to be compared
    error_for_comparison(optim_iter) = model.err;
end

 % find run that yieled smallest error value
[~, errorValsMin] = nanmin(error_for_comparison);

% save the run that gives the global min
modelData = monte(errorValsMin).model;

plot_vectorMod(modelData)
plot_vectorMod_data(modelData)
plot_vectorMod_model(modelData)


figure
pathPlot_hd(P, ST, MD); hold on;
plot(modelData.bestParams.xref, modelData.bestParams.yref, 'o')


