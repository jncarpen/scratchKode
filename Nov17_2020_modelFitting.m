% NOVEMBER 17, 2020: FIT COSINE MODEL

%% get cell info

% for simulated cell
for sessNum = 6:length(pos_cm)
    disp(strcat('processing session # ', sprintf('%.f', sessNum), '...'))
    if ~isempty(SpikeTimes{1,sessNum})
        for unitNum = 1:length(SpikeTimes{1,sessNum})
            disp(strcat('processing unit # ', sprintf('%.f', unitNum), '...'))
            if ~isempty(SpikeTimes{1,sessNum}{1,unitNum})
                P = pos_cm{1,sessNum};
                UID = UniqueID{1,sessNum}{1,unitNum};
                ST = SpikeTimes{1,sessNum}{1,unitNum};
                ref_point = hwCoord{1,sessNum};

                % angle_of_interest = 90;
                % [ST, ~, ~] = simulate_ego_cell(P, ref_point, angle_of_interest);

                %% fit model
                [model] = modelMe(P, ST);
                %cosFit or %cosErr


                %% plot 
                fig = figure('units','normalized','outerposition',[0 0 1 1]);
                set(gcf,'color','w');
                figtit = strcat('UID', sprintf('%.f', UID), '/S', sprintf('%.f', sessNum));
                sgtitle(figtit, 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')

                subplot(1,3,1)
                hold on;
                pathPlot_quiver(P, ST, get_hd(P));
                scatter(ref_point(1,1),ref_point(1,2),[70], 'b', 'filled')
                scatter(model.xref,model.yref,[70], 'r', 'filled')

                subplot(1,3,2)
                plot_modelfit(model)
        %         xline(deg2rad(angle_of_interest-180));

                subplot(1,3,3)
                plot_vectorMod(model)

                fileBody = strcat('UID', sprintf('%.f', UID), '_S', sprintf('%.f', sessNum));
                filename = strcat('D:\egoAnalysis\Nov17_vectorMod\', fileBody, '.png');
                saveas(fig, filename);
                
                close all
            end
        end
    end
end
