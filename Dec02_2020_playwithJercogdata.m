%% December 2, 2020
% play with the Jercog data set to figure out if I can reproduce any
% meaningful figures from the paper

%% LOAD DATA

% load dataset (.mat file)
load('D:\Data\KandelData\RatesNew10BinsEnvABC.mat')


%% FIGURE 2C
% read in excel sheet as an array
F2C = table2array(readtable('D:\Data\KandelData\FigureData.xlsx','Sheet','Figure 2c','Range','A2:B672'));

% variance explained by place tuning (distribution)
var_exp_place = F2C(:,1);
% variance explained by RH-angle model (distribution)
var_exp_model = F2C(:,2);

% plot
figure
subplot(1,2,1)
histogram(var_exp_place, 'FaceAlpha', .5, 'FaceColor', 'r');box off;
title('var explained by place', 'FontSize', 12); 
xlim([0 1]); xticks([0 .2 .4 .6 .8 1]); ylim([0 170])
subplot(1,2,2)
histogram(var_exp_model, 'FaceAlpha', .5, 'FaceColor', 'k'); box off;
title('var explained by model', 'FontSize', 12)
xlim([0 1]); xticks([0 .2 .4 .6 .8 1]); ylim([0 170])

%% FIGURE 3
F3 = table2array(readtable('D:\Data\KandelData\FigureData.xlsx','Sheet','Figure 3a,b,c,d,e,f','Range','A2:D698'));
g_vals = F3(:,1);
theta_vals = F3(:,2);
x_vals = F3(:,3); y_vals = F3(:,4);
max_bin = 10;

% 'pull in' the distant points
for d = 1:length(x_vals)
    if x_vals(d) < 0 % if negative
        x_vals(d) = -1;
    elseif x_vals(d) > 10
        x_vals(d) = 11;
    end
end
for d = 1:length(y_vals)
    if y_vals(d) < 0 % if negative
        y_vals(d) = -1;
    elseif y_vals(d) > 10
        y_vals(d) = 11;
    end
end


% plot
figure
subplot(2,2,1)
histogram(g_vals, 'FaceAlpha', .5, 'FaceColor', 'k');box off;
title('amplitude modulation', 'FontSize', 12)
subplot(2,2,2)
histogram(theta_vals, 'FaceAlpha', .5, 'FaceColor', 'k');box off;
title('RH-angle', 'FontSize', 12)
subplot(2,2,3)
histogram(x_vals, 'FaceAlpha', .5, 'FaceColor', 'k');box off;
title('x-reference', 'FontSize', 12); xlim([-2 12]); xticks([-1:2:11])
xticklabels({'distant', '1', '3', '5', '7', '9', 'distant'}); xtickangle(45);
subplot(2,2,4)
histogram(y_vals, 'FaceAlpha', .5, 'FaceColor', 'k');box off;
title('y-reference', 'FontSize', 12); xlim([-2 12]); xticks([-1:2:11])
xticklabels({'distant', '1', '3', '5', '7', '9', 'distant'}); xtickangle(45);


%%



