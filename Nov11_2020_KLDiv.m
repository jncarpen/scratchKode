for sessNum = 2%:length(SpikeTimes)
    
    % get position data for session
    position = pos_cm{1,sessNum};
    head_direction = hd{1,sessNum};
    x = position(:,2); y = position(:,3);
    t = position(:,1); fs = mode(diff(t));

    % choose reference points
    ref_point = hwCoord{1,sessNum};
    ref_point2 = boxCtr{1,sessNum};
    xGoal = ref_point(1,1);
    yGoal = ref_point(1,2);

    % get bin center locations
    nBins = 10; % divide the arena into 100 2D spatial bins
    [binCtrs] = get_spatial_bin_centers(position, nBins);

    for unitNum = 1:length(SpikeTimes{1,sessNum})
        S = SpikeTimes{1,sessNum}{1,unitNum};
        SpkTrn = SpikeTrain{1,sessNum}{1,unitNum};
        
        f = figure;
        set(gcf,'color','w');
        subplot(1,2,1)
        pathPlot_hd(position, S, head_direction)
        c1 = colorbar('hsv')
        cbfreeze(c1)
        h1 = plot(ref_point(1,1), ref_point(1,2), 'o', 'MarkerSize', 6);
        set(h1, 'markerfacecolor', 'k');
%         hold off;
        
        [hd_occ, allo_occ, ego_occ, time_occ, KL, KLI] = get_binned_occupancy(position, ref_point, "hd");
        subplot(1,2,2)
        imagesc(flipud(KL))
%         colormap(flipud(bone))
        cbfreeze('jet')
        colorbar
        c = colorbar;
        c.FontSize = 12;
        c.Box = "off";  c.FontName='Calibri';
        pbaspect([1 1 1])
        axis off
        title("Kullback-Leibler Divergence", 'FontSize', 14, 'FontName', 'Calibri')

    end
end