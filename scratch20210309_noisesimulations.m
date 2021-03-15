%% March 9, 2021 (Tues)
% Make simulted populations for noise calculations

%% 4. EGOCENTRIC BEARING
sesslist = [1:length(openfield), floor(linspace(4,84,14))];
for iter = 1:100
    % set parameters for hd component
    param = ego100(iter).param;
    param.noise = 2;
    param.P = openfield(sesslist(iter)).P;
    param.SNum = sesslist(iter); % save session #
    param.Z = get_hd(param.P);
    [sim, vm_pdf] = simulate_ego(param);
    ST = sim.ST;
    % pathPlot_hd(param.P, ST, param.Z);
    
    % run model
    modelData = modelMe(param.P, ST, param.Z);
    
    % save
    ego100_noise200(iter).param = param;
    ego100_noise200(iter).ST = ST;
    ego100_noise200(iter).out = modelData;
end
save('D:\Data\External Data\NOISE\noise-ego\ego100_noise200.mat', 'ego100_noise200','-v7.3');

