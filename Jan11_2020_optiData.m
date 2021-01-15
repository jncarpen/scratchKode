% January 11, 2020
% Start working with the OptiTrack data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data
% load('D:\Data\Dataset\25398\25398OT.mat');

for sess = 9:length(ST)
    % information for this session
    now = optiData{1,sess};
    
    if ~isempty(now)
        for unit = 1:length(ST{1,sess})
            if length(ST{1,sess}{1,unit}) > 100
                % grab the spikes now
                spikesNow = ST{1,sess}{1,unit};

                % yaw pitch and roll tuning curves
                [tc] = tc_YPR(now, spikesNow);

                %% plot
                titleNow = strcat('A25398', {' '}, 'S', sprintf('%.f', sess), ...
                    {' '}, 'U', sprintf('%.f', unit));
                fig = figure; set(gcf,'color','w');
                hold on;
                plot(tc.binCtrs, tc.yaw, 'r', 'LineWidth', 1.10)
                plot(tc.binCtrs, tc.pitch, 'k', 'LineWidth', 1.10)
                plot(tc.binCtrs, tc.roll, 'Color',[.6 .6 .6] , 'LineWidth', 1.10)
                xlim([-pi pi]); title(titleNow{1,1}, 'FontWeight', 'normal');
                ylabel("fr (Hz)"); xlabel("angle (rad)");
                l = legend('yaw', 'pitch', 'roll'); legend boxoff;
                set(gca,'FontSize',20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
                box off
                hold off;

                %% save stuff
                fileBody = strcat('S', sprintf('%.f', sess), ...
                '_U', sprintf('%.f', unit));

                filename = strcat('D:\egoAnalysis\Jan11_tcYPR\', fileBody, '.png');
                saveas(fig, filename);

                close all
            end
            
        end
        
    end
end








