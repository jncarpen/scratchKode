%% March 1, 2021
% Simulate example populations for the paper

%% LOAD TRAJECTORY DATA
% 100 open field sessions from jan sigurd's dataset
load('D:\Data\External Data\Blackstad-OF\dataof.mat')

%% 1. EGOCENTRIC BEARING + DISTANCE
for iter = 1:100
    param.A = (300-250).*rand(1) + 250;
    param.P = dataof(iter).P;
    param.Z = get_hd(param.P);
    param.theta = randi(360);
    param.kappa = 4;
    param.sigma = 8;
    param.radius = (35-25).*rand(1) + 25;
    param.rp = [(120-40).*rand(1) + 40,(120-40).*rand(1) + 40];
    % simulate cell & save spiketimes
    [sim] = simulate_egodist(param); ST = sim.ST;
    
    % monte carlo runs
    totalruns = 100;
    for run = 1:totalruns
        out = modelMe(param.P, ST, param.Z);
        monte(run).model = out; 
        err4compare(run) = out.model.error;
    end
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(err4compare);
    % save the run that gives the global min
    modelData = monte(errorValsMin).model;
    
    % save
    ring100(iter).param = param;
    ring100(iter).ST = ST;
    ring100(iter).out = modelData;
end
% save matfiles
save('D:\Data\External Data\Simulation-POP\ring100.mat', 'ring100', '-v7.3');

%% 2. PLACE 
for iter = 1:100
    root.A = (8-6).*rand(1) + 6;
    root.P = dataof(iter).P;
    root.Z = get_hd(root.P);
    root.sigma = [(15-10).*rand(1) + 10, (15-10).*rand(1) + 10];
    root.ctr = [(120-40).*rand(1) + 40,(120-40).*rand(1) + 40]; % exclude edges
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
    place100(iter).param = root;
    place100(iter).ST = ST;
    place100(iter).out = modelData;
end
save('D:\Data\External Data\Simulation-POP\place100.mat', 'place100', '-v7.3');


%% 3. PLACE + HD
for iter = 1:100
    root.A = (20-15).*rand(1) + 15;
    root.P = dataof(iter).P;
    root.Z = get_hd(root.P);
    root.sigma = [(15-12).*rand(1) + 12, (15-12).*rand(1) + 12];
    root.ctr = [(120-40).*rand(1) + 40,(120-40).*rand(1) + 40]; % exclude edges
    root.bins = 20;
    root.size = 150;
    
    % simulate cell & save spiketimes
    [map] = simulate_ratemap(root);
    % set parameters for hd component
    param.A = (30-25).*rand(1) + 25;
    param.P = root.P;
    param.Z = get_hd(root.P);
    param.kappa = (8-5).*rand(1) + 5;
    param.theta = randi(360);
    [sim] = simulate_placehd(map,param);
    ST = sim.ST;
    
    % monte carlo runs
    totalruns = 100;
    for run = 1:totalruns
        out = modelMe(param.P, ST, param.Z);
        monte(run).model = out; 
        err4compare(run) = out.model.error;
    end
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(err4compare);
    % save the run that gives the global min
    modelData = monte(errorValsMin).model;
    
    % save
    placehd100(iter).param = root;
    placehd100(iter).paramhd = param;
    placehd100(iter).ST = ST;
    placehd100(iter).out = modelData;
end
save('D:\Data\External Data\Simulation-POP\placehd100.mat', 'placehd100','-v7.3');


%% 4. EGOCENTRIC BEARING
for iter = 1:100
    % set parameters for hd component
    param.P = dataof(iter).P;
    param.A = (10-8).*rand(1) + 8;
    param.rp = [(120-40).*rand(1) + 40,(120-40).*rand(1) + 40]; 
    param.P = root.P;
    param.Z = get_hd(root.P);
    param.kappa = (8-5).*rand(1) + 5;
    param.theta = randi(360);
    [sim, vm_pdf] = simulate_ego(param);
    ST = sim.ST;
%     pathPlot_hd(param.P, ST, param.Z);
    
    % monte carlo runs
    totalruns = 100;
    for run = 1:totalruns
        out = modelMe(param.P, ST, param.Z);
        monte(run).model = out; 
        err4compare(run) = out.model.error;
    end
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(err4compare);
    % save the run that gives the global min
    modelData = monte(errorValsMin).model;
    
    % save
    ego100(iter).param = param;
    ego100(iter).ST = ST;
    ego100(iter).out = modelData;
end
save('D:\Data\External Data\Simulation-POP\ego100.mat', 'ego100','-v7.3');


%% 5. HEAD DIRECTION
for iter = 1:100
    % set parameters for hd component
    param.P = dataof(iter).P;
    param.A = (9-6).*rand(1) + 6;
    param.P = root.P;
    param.Z = get_hd(root.P);
    param.kappa = (8-5).*rand(1) + 5;
    param.theta = randi(360);
    [sim] = simulate_hd(param);
    ST = sim.ST;
%     pathPlot_hd(param.P, ST, param.Z);
    
    % monte carlo runs
    totalruns = 100;
    for run = 1:totalruns
        out = modelMe(param.P, ST, param.Z);
        monte(run).model = out; 
        err4compare(run) = out.model.error;
    end
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(err4compare);
    % save the run that gives the global min
    modelData = monte(errorValsMin).model;
    
    % save
    hd100(iter).param = param;
    hd100(iter).ST = ST;
    hd100(iter).out = modelData;
end
save('D:\Data\External Data\Simulation-POP\hd100.mat', 'hd100','-v7.3');

%% 6. PLACE + EGO (CTR)
for iter = 1:100
    root.A = (20-15).*rand(1) + 15;
    root.P = dataof(iter).P;
    root.Z = get_hd(root.P);
    root.sigma = [(15-12).*rand(1) + 12, (15-12).*rand(1) + 12];
    root.ctr = [(120-40).*rand(1) + 40,(120-40).*rand(1) + 40]; % exclude edges
    root.bins = 20;
    root.size = 150;
    
    % simulate cell & save spiketimes
    [map] = simulate_ratemap(root);
    % set parameters for hd component
    param.A = (70-60).*rand(1) + 60;
    param.P = root.P;
    param.Z = get_hd(root.P);
    param.kappa = 10; %(6-5).*rand(1) + 5;
    param.theta = randi(360);
    param.rp = root.ctr;
    [sim] = simulate_placeego(map,param);
    ST = sim.ST;
%     pathPlot_hd(param.P, ST, param.Z);
    
    % monte carlo runs
    totalruns = 100;
    for run = 1:totalruns
        out = modelMe(param.P, ST, param.Z);
        monte(run).model = out; 
        err4compare(run) = out.model.error;
    end
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(err4compare);
    % save the run that gives the global min
    modelData = monte(errorValsMin).model;
    
    % save
    placeego100(iter).param = root;
    placeego100(iter).paramhd = param;
    placeego100(iter).ST = ST;
    placeego100(iter).out = modelData;
end
save('D:\Data\External Data\Simulation-POP\placeego100.mat', 'placeego100','-v7.3');


