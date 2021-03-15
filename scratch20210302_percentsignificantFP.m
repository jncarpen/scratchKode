%% PERCENT SIGNIFICANT (FROM FALSE POSITIVE CALCULATIONS)
% March 2, 2021

%% pull out
clear all;
pathname = 'D:\Data\External Data\Trajectory-POP\FP\TJ63';
D = dir(pathname); D = D(3:102);

usortpathlist = cell(1,100);
for i = 1:100
    usortpathlist{i} = [D(i).folder, '/', D(i).name];
end
pathlist = natsortfiles(usortpathlist);

for i = 1:length(D)
    a = load(pathlist{i});
    tj63fp(i).B = a.B;
end

% save('D:\Data\External Data\NOISE\noise-ego\FP\ego100_noise0FP.mat', ...
%     'ego100_noise0FP', '-v7.3');

save('D:\Data\External Data\Trajectory-POP\FP\FPconcat\tj63fp.mat', 'tj63fp', '-v7.3')
load('D:\Data\External Data\Trajectory-POP\OFS63.mat')

%%

A = place100;
B = tj63fp;
%% find percentage of cells that are significant
HD_sig_NP = zeros(100,1); RH_sig_NP = zeros(100,1); 
for i = 1:100
    % grab the modulation strengths for each unit
    MS_HD_now = A(i).out.measures.TS.HD;
    MS_RH_now = A(i).out.measures.TS.RH;
    
    % grab the shuffled distribution for each unit
    HD_shuf_now = B(i).B.mshd;
    RH_shuf_now = B(i).B.msrh;
    
    % find confidence interval for the shuffled distribution
    HD_ci = [mean(HD_shuf_now) - 2*std(HD_shuf_now), mean(HD_shuf_now) + 2*std(HD_shuf_now)];
    RH_ci = [mean(RH_shuf_now) - 2*std(RH_shuf_now), mean(RH_shuf_now) + 2*std(RH_shuf_now)];
    
    HD_sig(i) = MS_HD_now < HD_ci(1) | MS_HD_now > HD_ci(2);
    RH_sig(i) = MS_RH_now < RH_ci(1) | MS_RH_now > RH_ci(2);
    
    % sort the shuffled distribution and real values
    HD_sort = sort([HD_shuf_now, MS_HD_now], 'ascend');
    RH_sort = sort([RH_shuf_now, MS_RH_now], 'ascend');
    
    % use a non-parametric method to determine significance
    if find(HD_sort == MS_HD_now) > 950 | find(HD_sort == MS_HD_now) < 50; HD_sig_NP(i) = 1; end
    if find(RH_sort == MS_RH_now) > 950 | find(RH_sort == MS_RH_now) < 50; RH_sig_NP(i) = 1; end
end

disp('COMPUTING MODEL SENSITIVITY MEASURES...')
disp(['HD TS (parametric):', num2str(sum(HD_sig)), '%']);
disp(['HD TS (non-parametric):', num2str(sum(HD_sig_NP)), '%']);
disp(['RH TS (parametric):', num2str(sum(RH_sig)), '%']);
disp(['RH TS (non-parametric):', num2str(sum(RH_sig_NP)), '%']);


%%


%% find values for each run
for i=1:100
    MS_HD(i) = A(i).out.measures.TS.HD;
    MS_RH(i) = A(i).out.measures.TS.RH;
    
    MS_HD_shuff(i) = mean(B(i).B.mshd, 'omitnan');
    MS_RH_shuff(i) = mean(B(i).B.msrh, 'omitnan');
    
    VE_place_shuff(i) = mean(B(i).B.vep, 'omitnan');
    VE_RH_shuff(i) = mean(B(i).B.verh, 'omitnan');
    
    VE_place(i) = A(i).out.measures.VE.place;
    VE_RH(i) = A(i).out.measures.VE.RH;
    
end


%% PLOT
plotName = 'RH VE';
X = MS_RH;
Xshuf = MS_RH_shuff;
stdev = std(X, 'omitnan');
figure;
set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', [0 .3 .8]);
histogram(Xshuf, 15, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'k');
xline(mean(Xshuf) + stdev*2, ':k', 'LineWidth', 1.5); % std
xline(mean(Xshuf) - stdev*2, ':k', 'LineWidth', 1.5); % std
ylabel("frequency (count)"); xlabel("tuning strength");
title(plotName, 'FontName', 'Helvetica UI', 'FontSize', 20, 'FontWeight', 'normal');
set(gca,'FontSize',20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
l = legend('real data', 'shuffled data', '95% CI', 'Location', 'northeastoutside');
legend boxoff    
box off;


