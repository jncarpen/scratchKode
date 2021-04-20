%% 2. PLACE 
clear place100
sess = 63;
for iter = 1:100
    root = A.place100(iter).param;
    root.P = dataof(sess).P;
    root.Z = get_hd(root.P);
    
    % simulate cell & save spiketimes
    [map] = simulate_ratemap(root);
    [sim] = simulate_place(map,root.P(:,1)); ST = sim.ST;
%     pathPlot_hd(root.P, ST, get_hd(root.P))
    
    modelData = modelMe(root.P, ST, get_hd(root.P));
    
    % save
    place100(iter).param = root;
    place100(iter).ST = ST;
    place100(iter).out = modelData;
end
save('D:\Data\External Data\Trajectory-POP\OFS63.mat', 'place100', '-v7.3');
