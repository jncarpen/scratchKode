
% clear all; clc;
% load in some data
data = load('D:\Data\Dataset\24116\24116_2.mat');

ofsess = [];
count = 1;
for i = 1:numel(data.pos_cm)
    now = data.trialType{1,i};
    if contains(now, 'OF')
        ofsess(count) = i;
        count = count + 1;
    end
end

% loop through all open field sessions
for j = 1:numel(ofsess)
        sess = ofsess(j);
        P = data.pos_cm{1,sess};
    for unit = 1:length(data.SpikeTimes{1,sess})
        ST = data.SpikeTimes{1,sess}{1,unit};
        figure;
        pathPlot_hd(P, ST, get_hd(P));
        title(['Sess:', num2str(sess), 'Unit: ', num2str(unit)]);
        pause; close all;
    end
end


% use the first session
sess = 14;
P = data.pos_cm{1,14};
param.theta = 180;
param.Z = get_hd(P)-180;
param.P = P;
param.kappa = 5;
param.rp = [75, 75];
param.A = 10;

[sim] = simulate_ego(param);
ST = sim.ST;

pathPlot_hd(P, ST, get_hd(P))


















