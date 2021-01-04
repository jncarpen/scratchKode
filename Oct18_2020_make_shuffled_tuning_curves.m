%% October 22, 2020

% for animal #25398
[SpikeTimes_thresh] = speedThreshold_loop(pos_cm, speed_cm, SpikeTimes);

for sessNum = 1:length(SpikeTrain)
     if ~isempty(SpikeTimes{1,sessNum}) % skip empty trials
        disp(strcat("processing session ",  sprintf('%.f', sessNum), "..."))
        
        % get information about current session
        P = pos_cm{1,sessNum};
        
        % get reference points
        ctrCoord = boxCtr{1,sessNum}; % box center
        refCoord = hwCoord{1,sessNum}; % home well
        
        for unit = 1:length(SpikeTimes{1,sessNum})
            if length(SpikeTimes_thresh{1,sessNum}{1,unit}) > 100
                disp(strcat("processing unit ",  sprintf('%.f', unit), "..."))
                
                % get spike times (thresholded)
                ST_thresh = SpikeTimes_thresh{1,sessNum}{unit};
                
                % get spike times (unthresholded)
                ST = SpikeTimes{1,sessNum}{unit};
                
                % get unique ID #
                UID = UniqueID{1,sessNum}{1,unit};
            
                %% plot
                fig = figure('units','normalized','outerposition',[.2 .2 .6 .65]);
                set(gcf,'color','w');
                
                figTit = strcat('UID', sprintf('%.f', UID), 'SESS', sprintf('%.f', sessNum));
                fig.Name = figTit; % set figure name
                sgtitle(figTit);
                
                % plot tuning curve
                subplot(1,2,1)
                tc_shuffle(P, ST_thresh, refCoord);
                
                % plot HD pathplot
                subplot(1,2,2)
                HD = get_hd(P);
                pathPlot_hd(P, ST_thresh, HD);
                
                % save figure
                filename = strcat('D:\egoAnalysis\Oct25_tc_shuffled_25398\', figTit, '.png');
                saveas(fig, filename);
                
                % get ready for next figure
                close all
            else
                % display message if there are not enough spikes to include
                message = strcat("unit ", sprintf('%.0f', unit), " has too few spikes to include");
                disp(message)
            end
     end
  end
end
