
clear  whosJZ_n  whosJZ_u  whosJZ_count

% count = 2206;

for n = 861:numel(JZM.neurons)
    for u = 1:numel(JZM.neurons(n).members)
        
%         whosJZ_n(count) = n;
%         whozJZ_u(count) = u;
%         whosJZ_count(count) = count;
%         count = count + 1;
        
        uid = JZ.neurons(n).members(u).uid;
        P_raw = units2cm(JZ.neurons(n).members(u).P);
        ST = JZ.neurons(n).members(u).ST;
        session = JZ.neurons(n).members(u).sessName;
        
        % smooth position
        clear P
        sigma = 2;
        P = smooth_pos(P_raw, sigma);
        
        model = JZM.neurons(n).members(u).model;
        ref_point = [model.bestParams.xref, model.bestParams.yref];
        
        ve_place = model.varExplained.place;
        ve_model = model.varExplained.model;
        g = model.bestParams.g;
        RH_angle = model.bestParams.thetaP;
        xx = model.bestParams.xref;
        yy = model.bestParams.yref;
        
        % plot
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color','w');
        fileBody = strcat('N:', sprintf('%.f', n), ', U:', sprintf('%.f', u), ', IX:', sprintf('%.f', count), ', VEP:', sprintf('%.2f', ve_place), ', VEM:', sprintf('%.2f', ve_model), ', g:', sprintf('%.f', g), ', theta p:', sprintf('%.f', RH_angle), ', X:', sprintf('%.1f', xx), ', Y:', sprintf('%.1f', yy));
        sgtitle(fileBody, 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
        
        
        if xx > 0 && xx < 150 && yy > 0 && yy < 150
            subplot(2,2,1)
            pathPlot_hd(P, ST, get_hd(P))
            hold on;
            scatter(ref_point(1,1), ref_point(1,2), 60, 'k', 'filled')            
        else
            subplot(2,2,1)
            pathPlot_hd(P, ST, get_hd(P))
        end
        hold off;
        
        
        subplot(2,2,2)
        [tcVals_egoAng, binCtrs_egoAng] = plot_egoBearing(P, ST, ref_point, "True")
        get(gca)
        set(gca,'FontSize',15, 'FontName', 'Calibri Light')
        pbaspect([1 1 1])
        
        if sum(sum(~isnan(model.rateMap))) > 2
            subplot(2,2,3)
            plot_vectorMod(model)
        end
        
        subplot(2,2,4)
        map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
        peakRate = nanmax(nanmax(map.z));
        rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
        plot.colorMap(map.z)
        pbaspect([1 1 1])
        colormap(gca,'jet')
        c2 = colorbar; c2.FontSize = 15; c2.FontName = 'Calibri Light'; c2.FontWeight = 'normal';
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 15, 'FontWeight', 'normal');
        box off
        
        fileBody = strcat('N', sprintf('%.f', n), '_U', sprintf('%.f', u), '_IX', sprintf('%.f', count));
        filename = strcat('D:\egoAnalysis\Dec11_dManConversion\', fileBody, '.png');
        saveas(fig, filename);
        close all
        
        count = count + 1;
    end
end



%% pull out useful information

count = 1;
for n = 1:numel(JZM.neurons)
    for u = 1:numel(JZM.neurons(n).members)
        
        modnow = JZM.neurons(n).members(u).model;
        gparam(count) = modnow.bestParams.g;
        varex(count) = modnow.varExplained.model;
        xparam(count) = modnow.bestParams.xref;
        yparam(count) = modnow.bestParams.yref;
        thetaparam(count) = modnow.bestParams.thetaP;
        err(count) = modnow.err;
        
        count = count + 1;
        
    end
end


% make histograms
% remove outliers
ve = thetaparam;

% make
figure
hold on;
myBins = 0:30:360
histogram(ve, myBins,'FaceColor', 'r');   
ax = gca; alpha(ax,.1);
title("Preferred Theta Parameters", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')



xlabel("angle (deg)", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')
ylabel("frequency", "FontSize", 20, 'FontName', "Calibri Light", 'FontWeight', 'bold')



x_vals = xparam; y_vals = yparam;
% 'pull in' the distant points
for d = 1:length(x_vals)
    if x_vals(d) < 0 % if negative
        x_vals(d) = -1;
    elseif x_vals(d) > 150
        x_vals(d) = 155;
    end
end
for d = 1:length(y_vals)
    if y_vals(d) < 0 % if negative
        y_vals(d) = -1;
    elseif y_vals(d) > 150
        y_vals(d) = 155;
    end
end

bins = -10:10:160
histogram(y_vals, bins, 'FaceAlpha', .5, 'FaceColor', 'k');box off;
title('x-reference', 'FontSize', 12); xlim([-2 155]); xticks([-10:50:170])
% xticklabels({'distant', '1', '3', '5', '7', '9', 'distant'}); xtickangle(45);


% 
figure
d = -pi:0.01:pi
aa = cos(d);
plot(d,aa, 'LineWidth', 1.5, 'Color', 'k'); hold on;
aa2 = -cos(d);
plot(d,aa2, '--k','LineWidth', 1.5);
bb = .5.*cos(d);
plot(d,bb, 'r', 'LineWidth', 1.5); 
cc = -.5.*cos(d);
plot(d,cc, '--r', 'LineWidth', 1.5); 
dd = 2.*cos(d);
plot(d,dd, 'b', 'LineWidth', 1.5); 
ee = -2.*cos(d)
plot(d,ee, '--b', 'LineWidth', 1.5); 
legend('cos(x)', '-cos(x)', '1/2cos(x)', '-1/2cos(x)', '2cos(x)', '-2cos(x)')
xlim([-pi pi])
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2', '\pi'})
box off;


pathPlot_hd(P, ST, get_hd(P))












%% SCRATCH
%         [boxCtrX,boxCtrY] = getBoxCenter(P);
%         ref_point = [boxCtrX,boxCtrY];
%         JZM.neurons(n).members(u).model = model;


%         
%         subplot(1,3,1)
%         plot_modelDynamics(P, ST, model, ref_point)
% 
%         subplot(1,3,2)
%         plot_iterations(model)
% 
%         subplot(1,3,3)
%         plot_vectorMod(model)












