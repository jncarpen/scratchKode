%%  simulate a bursty neuron
%   j. carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% randomly choose a session
% whichSession = randi(length(JZ.neurons));
% whichUnit = length(JZ.neurons(whichSession).members);

whichSession = 1232;
whichUnit = 2;

% grab data from that session
P_raw = units2cm(JZ.neurons(whichSession).members(1).P);
P = smooth_pos(P_raw, 2);
ST = JZ.neurons(whichSession).members.ST;

% plot to visualize real data
figure(1)
pathPlot_hd(P, ST, get_hd(P))

% ratemap
map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 3);
figure; imagesc(map.z);

% simulate bursty cell
[sim] = simulate_bursty_place(map, P);

% run the model
initial = choose_initial_conditions(P);
[model] = modelMe(P, sim.ST, get_hd(P), initial);


% visualize simulated cell
figure(2)
pathPlot_hd(P, sim.ST, get_hd(P))

figure(3)
map = analyses.map(P, sim.ST, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
colormap(gca,'jet')
c2 = colorbar; c2.FontSize = 25;
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
box off;

% loop through a bunch of neurons
for i = 1:length(neuron)
    now = neuron(i);
    ST_now = JZ.neurons(now).members(1).ST;
    
    if length(ST_now) > 100
        % grab position information now; smooth & convert to cm
        P_now = smooth_pos(units2cm(JZ.neurons(now).members(1).P),2);
        
        % make the ratemap
        map = analyses.map(P_now, ST_now, 'smooth', 2, 'binWidth', 3);
        
        % simulate a bursty neuron
        [sim] = simulate_bursty_place(map, P_now);
        
        % run the model
        initial = choose_initial_conditions(P_now);
        [model] = modelMe(P_now, sim.ST, get_hd(P_now), initial);
    end
    
end
    






