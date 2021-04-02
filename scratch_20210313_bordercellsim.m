%% MAKE A POPULATION OF PLACE CELLS (BORDER CELLS)
% Larger, skewed fields with place field centers on the edges of box
tic
addpath(genpath('/cluster/home/jordannc/data'))
addpath(genpath('/cluster/home/jordannc/data/simulate_spatial'))

OF = load('/cluster/home/jordannc/data/Simulation-POP/openfield.mat');
OF = OF.openfield;
sesslist = [1:length(OF), floor(linspace(4,84,14))];
ctrlist = [40*rand(100,1);(150-120).*rand(100,1) + 120];

for iter = 1:100
    root.P = OF(sesslist(iter)).P;
    root.Z = get_hd(root.P);
    root.A = (8-6).*rand(1) + 6;
    root.sigma = [(18-12).*rand(1) + 12, (18-12).*rand(1) + 12];
    root.ctr = [datasample(ctrlist,1),datasample(ctrlist,1)];
    root.bins = 20;
    root.size = 150;
    
    % simulate cell & save spiketimes
    [map] = simulate_ratemap(root);
    [sim] = simulate_place(map,root.P(:,1)); ST = sim.ST;
    pathPlot_hd(root.P, ST, get_hd(root.P))
        
    % monte carlo runs
    totalruns = 100;
    for run = 1:totalruns
        out = modelMe(root.P, ST, get_hd(root.P));
        monte(run).model = out; 
        err4compare(run) = out.model.error;
    end
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(err4compare);
    % save the run that gives the global min
    modelData = monte(errorValsMin).model;
    
    % save
    border100(iter).param = root;
    border100(iter).ST = ST;
    border100(iter).out = modelData;
end

save('D:\Data\External Data\Trajectory-POP\OF_S2\border100.mat', 'border100', '-v7.3');
