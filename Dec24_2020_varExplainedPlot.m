%% December 24, 2020
% Check to see if all variance explained values are between 0 and 1 and if 
% they can be plotted as they are in the Jercog paper

clear varex_place varex_model modstren_hd modstren_rh err test
count = 1;
for nn = 752:1252 %length(JZ.neurons)
    for uu = 1:length(JZ.neurons(nn).members)
        
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(nn).members(uu).P);
        ST = JZ.neurons(nn).members(uu).ST;

        % smooth position vectors
        clear P
        sigma = 2;
        P = smooth_pos(P_raw, sigma);
        HD = get_hd(P);
        
        initial = choose_initial_conditions(P);
        [model] = modelMe(P, ST, HD, initial);
        
        varex_place(count) = model.varExplained.place;
        varex_model(count) = model.varExplained.model;
        
        modstren_hd(count) = model.modStrength.HD;
        modstren_rh(count) = model.modStrength.RH;
        
        err(count) = model.err;
        
        test(count).model = model; % save the model for each iter
        test(count).neuron = nn;
        test(count).unit = uu;
        test(count).ST = ST;
        test(count).P = P;
        test(count).HD = HD;

        count = count + 1;

    end
end


% figure out why some of the variance explained values are negative
for i = 1:length(err)
    if varex_model(i) < 0 
%         model = test(i).model;
%         ref_point = [model.bestParams.xref,model.bestParams.yref];
        P = test(i).P; ST = test(i).ST;
        pathPlot_hd(P, ST, get_hd(P))
        pause
        close all
        
%         t = P(:,1);
        
%         length(ST)

%         % get angles associated with each spike
%         spkhd = HD(knnsearch(t, ST));
% 
%         % make tuning curve (smooth and big bins)
%         tc = analyses.turningCurve(spkhd, HD, Fs, 'smooth', 3, 'binWidth', 3);  
%         
%         % plot  
%         fig = figure;
%         plot(tc(:,1), tc(:,2))
%         
%         fileBody = strcat('unit', sprintf('%.f', i));
%         filename = strcat('D:\egoAnalysis\Dec24_negVarEx_HDtc\', fileBody, '.png');
%         saveas(fig, filename);
%         
%         close all
    end
    
end




% name the variables
X = varex_model;
plotName = 'variance explained by model';
xname = 'variance explained';

% plot
figure; set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'FaceColor', [1 .2 0]);
ylabel("frequency (count)"); xlabel(xname);
title(plotName);
set(gca,'FontSize',20, 'FontName', 'Calibri Light', 'FontWeight', 'normal');
box off;

