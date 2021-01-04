%% DECEMBER 14, 2020
% Make figures for Jan Sigurd's Lab Meeting Presentation (Dec 15, 2020)

% manually define neurons/units of interest
% n_interest = [524 495 467 468 459 457 440 285 284 283 282];
% u_interest = [6 3 2 1 1 1 1 2 2 2 4];

n_interest = [1196 1140 1077 1079 1064 1060 1035 798 794 792 788 597 587 560 542 498 486 420];
% decide on how many mc simulations to run
total_iters = 1;


clear errorValsMin errorVals dclub2
count_neuron = 1;
for nn = 1:length(n_interest)
    n_now = n_interest(nn);
    for uu = 1:length(JZ.neurons(n_now).members)
        
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(n_now).members(uu).P);
        ST = JZ.neurons(n_now).members(uu).ST;

        % smooth position vectors
        clear P
        sigma = 2;
        P = smooth_pos(P_raw, sigma);
        HD = get_hd(P);

        % monte carlo
        for ii = 1:total_iters
            % generate random values for parameter intial conditions
            initial = choose_initial_conditions(P);
            % run the model
            [model] = modelMe(P, ST, HD, initial);
%             dclub.neurons(cellcell).model(ii) = model;
%             dclub.neurons(cellcell).initial(ii) = initial;
            
            dclub2.neurons(nn).members(uu).model(ii) = model;
            dclub2.neurons(nn).members(uu).initial(ii) = initial;
            
            errorVal_unit(ii) = model.err;
        end

        % find the tiniest error value across the 100 simulations
        [~, errorValsMin{1,count_neuron}(uu)] = nanmin(errorVal_unit);
    end
    count_neuron = count_neuron + 1;
end











% now that we have found the iteration that retrieved the 
% global minimum ('errorValsMin') we can plot the results -->
%% FIGURES

count = 1;
for nn = 1:length(n_interest)
    n_now = n_interest(nn);
    
    for uu = 1:length(JZ.neurons(n_now).members)
        % which neuron/unit is being processed
        disp(strcat("nn = ", sprintf('%.f', nn), " and uu = ", sprintf('%.f', uu)))        
        
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(n_now).members(uu).P);
        ST = JZ.neurons(n_now).members(uu).ST;

        % smooth position
        clear P
        sigma = 2;
        P = smooth_pos(P_raw, sigma);

        % grab model that yielded global min
        tinyError = errorValsMin{1,nn}(uu);
        model = dclub2.neurons(nn).members(uu).model(tinyError);
        initVals = dclub2.neurons(nn).members(uu).initial(tinyError);
        % model = dclub.neurons(cellcell).model(tinyError);
        % initVals = dclub.neurons(cellcell).initial(tinyError);

        % find reference point fitted by model
        ref_point = [model.bestParams.xref, model.bestParams.yref];
        
        % get information about each neuron
%         info_stuff(count, 1) = nn;
%         info_stuff(count, 2) = n_now;
%         info_stuff(count,3) = uu;
%         info_stuff(count, 4) = model.bestParams.g;
%         info_stuff(count, 5) = model.bestParams.thetaP;
%         info_stuff(count, 6) = model.bestParams.xref;
%         info_stuff(count, 7) = model.bestParams.yref;
%         

%         if ref_point(1,1) > -10 && ref_point(1,1) < 180 && ref_point(1,2) > -10 && ref_point(1,2) < 180
            if sum(sum(~isnan(model.rateMap))) > 3
                fig = figure; set(gcf,'color','w');
                subplot(1,1,1);
                plot_vectorMod(model)

            else
                fig = figure; set(gcf,'color','w');
                subplot(1,1,1);
            end
                
            
%             fig = figure('units','normalized','outerposition',[0 0 1 1]);
%             set(gcf,'color','w');
%             subplot(1,1,1)
% %             tc_shuffle(P, ST, ref_point)
%             
%             egoRateMap(P, ST, ref_point)
            
%             pathPlot_hd(P, ST, get_hd(P)); hold on;
%             plot(ref_point(1,1), ref_point(1,2), 'o', 'MarkerFaceColor', 'k', 'MarkerSize', [13])
            
%         else
%             fig = figure; set(gcf,'color','w');
%             subplot(1,1,1)
%             pathPlot_hd(P, ST, get_hd(P)); hold on;
            
%         end


%         save
        fileBody = strcat('n_', sprintf('%.f', n_now), 'u_', sprintf('%.f', uu));
        filename = strcat('D:\egoAnalysis\Dec14_JS_labmeeting2\vectorMod2\', fileBody, '.png');
        saveas(fig, filename);

        close all
        
    end

end


% make ratemap

figure
map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
colormap(gca,'jet')
c2 = colorbar; c2.FontSize = 25;
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
box off









