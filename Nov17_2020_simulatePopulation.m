% November 17, 2020
% make a dataset of simulated cells

clear position refPoints ST spikes

position = pos_cm;
refPoints = hwCoord;
ST_sim = cell(1, length(pos_cm));

for sessNum = 1:length(pos_cm)
    P = pos_cm{1,sessNum};
    ref_point = hwCoord{1,sessNum};
    clear spikes
    spikes = cell(1, 8);
    count = 1;
    for AOI = [45:45:360]
        [ST, ~, ~] = simulate_ego_cell(P, ref_point, AOI);
        spikes{1,count} = ST;
        count = count + 1;
    end
    ST_sim{1,sessNum} = spikes;
end

% save in a struct
sim.position = position;
sim.refPoints = refPoints;
sim.ST = ST_sim;


%% Run the model on the simulated dataset
angle_of_interest = [45:45:360];
for sessNum = 50%:length(sim.position)
    for unitNum = 1%:length(sim.ST{1,sessNum})
        % grab information for session
        P = sim.position{1,sessNum};
        ST = sim.ST{1,sessNum}{1,unitNum};
        ref_point = sim.refPoints{1,sessNum};
        AOI = angle_of_interest(unitNum);
        
        % fit the model
        [model] = modelMe(P, ST)
        %% plot 
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color','w');
        figtit = strcat('S ', sprintf('%.f', sessNum), ', AOI ', sprintf('%.f', AOI), '(deg)');
        sgtitle(figtit, 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')

        subplot(1,3,1)
        hold on;
        pathPlot_quiver(P, ST, get_hd(P));
        scatter(ref_point(1,1),ref_point(1,2),[70], 'b', 'filled')
        scatter(model.xref,model.yref,[70], 'r', 'filled')

        subplot(1,3,2)
        plot_modelfit(model)
        xline(deg2rad(AOI-180));

        subplot(1,3,3)
        plot_vectorMod(model)

        fileBody = strcat('S', sprintf('%.f', sessNum), '_AOI', sprintf('%.f', AOI));
        filename = strcat('D:\egoAnalysis\Nov17_vectorMod_sim\', fileBody, '.png');
        saveas(fig, filename);

        close all
        
        
    end
end




