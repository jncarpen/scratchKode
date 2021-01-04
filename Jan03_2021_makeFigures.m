% make figures

% load dataset from the first animal
load("D:\Data\Dataset\24116\24116_2.mat")

% choose position information from a random session
P = pos_cm{1, randi(length(pos_cm))};
possible_x = min(P(:,2))+10:0.5:max(P(:,2))-10;
possible_y = min(P(:,3))+10:0.5:max(P(:,3))-10;

% set parameters for simulation
clear param
param.position = P;
param.ctr_mass = [rand(1)*150, rand(1)*150];
param.noise = 0;
param.width = (rand(1)*15)+6;

% simulate a place cell
[sim] = simulate_place(param);

% handle outputs
ST = sim.spiketimes;
spkTrn = binSpikes(P(:,1), ST);
HD = get_hd(P);

% perform the optimization (using a monte carlo method)
total_iters = 100; 
clear monte error_for_comparison
for optim_iter = 1:total_iters 
    % choose some initial conditions randomly
    initial = choose_initial_conditions(P);
    % run the model 
    [model] = modelMe(P, ST, HD, initial);

    monte(optim_iter).model = model; % save all the runs to be compared
    error_for_comparison(optim_iter) = model.err;
end

% find run that yieled smallest error value
[~, errorValsMin] = nanmin(error_for_comparison);

% save the run that gives the global min
bestModel = monte(errorValsMin).model;

% ref_point as chosen by model
ref_point = [bestModel.bestParams.xref, bestModel.bestParams.yref];

% plot stuff
plot_vectorMod(bestModel)

% make a 'root' struct for egocentric rate map
clear root
root.x = P(:,2);
root.y = P(:,3);
root.md = deg2rad(HD); % in radians
root.ts = P(:,1);
root.spike = spkTrn;

% run egocentric rate map
[out] = EgocentricRatemap(root, 'distanceBins', 0:1:100);
fign = 1;
plotEBC(root,out,fign)


