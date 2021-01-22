for i = 1:length(pos_cm)
    P_now = pos_cm{1,i};
    EEG_now = rawEEG{1,i}{1,1};
    
    for unit = 1:length(SpikeTimes{1,i})
        ST_now = SpikeTimes{1,i}{1,unit};
        % plot
        figure
        title_now = strcat('N', sprintf('%.f', i), ' U', sprintf('%.f', unit));
        pathPlot_theta(P_now, ST_now, EEG_now);
        title(title_now);
        pause; close all;
    end
end

pathPlot_theta(pos_, SpikeTimes_, rawEEG)