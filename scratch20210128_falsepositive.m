%% FALSE POSITIVE RATE (WITH SIMULATED CELLS)
%   January 25, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% load Pablo's data
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example1_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example2_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example3_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example4_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example5_RawData.mat')

% arrange dataset
Jercog(1).P = Cell1_Trajectory;
Jercog(2).P = Cell2_Trajectory;
Jercog(3).P = Cell3_Trajectory;
Jercog(4).P = Cell4_Trajectory;
Jercog(5).P = Cell5_Trajectory;

warning('off','all')
clear Msim
for simcell = 1:100
    % tell user which iteration youre on
    dispMe = strcat('running simulation ', {' '}, sprintf('%.f', simcell), ' of 100...');
    disp(dispMe{1,1});
    
    % grab some position data
    whichSession = randi(5);
    P = Jercog(whichSession).P;
    P(:,2:end) = P(:,2:end)+60;
    Z = rad2deg(get_MD(P))';
    
    % set maximum arena size (cm)
    maxSz = 120;
    minSz = 0;
    
    % set ratemap parameters
    root.peakrate = 10;
    root.ctr = [((maxSz-maxSz*.1)-(minSz+maxSz*.1)).*rand(1) + minSz+maxSz*.1,...
        ((maxSz-maxSz*.1)-(minSz+maxSz*.1)).*rand(1) + minSz+maxSz*.1];
    root.sigma = [(20-15).*rand(1) + 15, (20-15).*rand(1) + 15];
    root.size = maxSz;
    root.bins = 50;
    root.P = P;
    
    % simulate place cell
    [map] = simulate_ratemap(root);
    [sim] = simulate_place(map, P(:,1));
    
    % spiketimes
    ST = sim.ST;
    
    % save some stuff
    Msim(simcell).ST = ST;
    Msim(simcell).param = root;

    % sampling frequency (50 samples/sec)
    Fs = mode(diff(P(:,1))); 

    % perform the optimization (using a monte carlo method)
    total_iters = 100; 
    clear monte error_for_comparison
    for optim_iter = 1:total_iters 
        % run the model 
        out = modelMe(P, ST, Z);
        % save all the runs to be compared
        monte(optim_iter).model = out; 
        error_for_comparison(optim_iter) = out.model.error;
    end

    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(error_for_comparison);

    % save the run that gives the global min
    modelData = monte(errorValsMin).model;

    % save stuff
    Msim(simcell).modelData = modelData;


    %% shuffle the head direction values
    % define the number of shuffles we want
    total_shuffles = 1000;
    min_shuffle = floor(60/Fs); % 30 seconds
    max_shuffle = floor(600/Fs); % 200 seconds

    % systematic shuffle (can be replaced with a pseudorandom shuffle)
    shuffle_vec = [linspace(-min_shuffle, -max_shuffle, total_shuffles/2), ...
        linspace(min_shuffle, max_shuffle, total_shuffles/2)]; 
    shuffle_vec =  floor(shuffle_vec); % make them all integers

    % clear stuff from before
    clear param_g param_thetaP param_xref param_yref ...
        varex_place varex_rh modstren_hd modstren_hd error_shuff...
        varex_place_sig varex_rh_sig modstren_rh_sig  modstren_hd_sig...
        error_sig

    % shuffle the head direction values 1000x
    for shuff_num = 1:total_shuffles
        clear model_shuffled
        % grab the current shuffle value from vector
        shuff_value = shuffle_vec(shuff_num);

        % circularly shifts the elements in array A by K positions
        shuffled_Z = circshift(Z, shuff_value);

        % run the model on the shuffled data
        [out_SH] = modelMe(P, ST, shuffled_Z);

        % SAVE OUTPUT FROM EACH RUN
        
        % (1) best fit parameters
        param_g(shuff_num) = out_SH.model.fitParams.g;
        param_thetaP(shuff_num) = out_SH.model.fitParams.thetaP;
        param_xref(shuff_num) = out_SH.model.fitParams.xref;
        param_yref(shuff_num) = out_SH.model.fitParams.yref;
        
        % (2) variance explained (place & model)
        ve_place(shuff_num) = out_SH.measures.VE.place;
        ve_rh(shuff_num) = out_SH.measures.VE.RH;
        
        % (3) modulation/tuning strength
        TS_hd(shuff_num) = out_SH.measures.TS.HD;
        TS_rh(shuff_num) = out_SH.measures.TS.RH;
        
        % (4) error
        error_shuff(shuff_num) = out_SH.model.error;

    end
    
    % save stuff
    Msim(simcell).ve_place = ve_place;
    Msim(simcell).ve_rh = ve_rh;
    Msim(simcell).ts_hd = TS_hd;
    Msim(simcell).ts_rh = TS_rh;
    Msim(simcell).shuffErr = error_shuff;
    Msim(simcell).g = param_g;
    Msim(simcell).thetap = param_thetaP;
    Msim(simcell).xref = param_xref;
    Msim(simcell).yref = param_yref;
end
warning('on','all')

%% find percentage of cells that are significant
% @todo: implement non-parametric method
% clear HD_sig RH_sig HD_sig_NP RH_sig_NP 
HD_sig_NP = zeros(100,1); RH_sig_NP = zeros(100,1); 
for i = 1:100
    % grab the modulation strengths for each unit
    MS_HD_now = Msim(i).modelData.measures.TS.HD;
    MS_RH_now = Msim(i).modelData.measures.TS.RH;
    
    % grab the shuffled distribution for each unit
    HD_shuf_now = Msim(i).ts_hd;
    RH_shuf_now = Msim(i).ts_rh;
    
    % find confidence interval for the shuffled distribution
    HD_ci = [mean(HD_shuf_now) - 2*std(HD_shuf_now), mean(HD_shuf_now) + 2*std(HD_shuf_now)];
    RH_ci = [mean(RH_shuf_now) - 2*std(RH_shuf_now), mean(RH_shuf_now) + 2*std(RH_shuf_now)];
    
    HD_sig(i) = MS_HD_now < HD_ci(1) | MS_HD_now > HD_ci(2);
    RH_sig(i) = MS_RH_now < RH_ci(1) | MS_RH_now > RH_ci(2);
    
    % sort the shuffled distribution and real values
    HD_sort = sort([HD_shuf_now, MS_HD_now], 'ascend');
    RH_sort = sort([RH_shuf_now, MS_RH_now], 'ascend');
    
    % use a non-parametric method to determine significance
    if find(HD_sort == MS_HD_now) > 950; HD_sig_NP(i) = 1; end
    if find(RH_sort == MS_RH_now) > 950; RH_sig_NP(i) = 1; end
end

sum(HD_sig)
sum(HD_sig_NP)

sum(RH_sig)
sum(RH_sig_NP)

%% find values for each run
clear MS_HD MS_RH MS_HD_shuff MS_RH_shuff ...
    VE_place VE_RH
for i=1:100
    MS_HD(i) = Msim(i).modelData.measures.TS.HD;
    MS_RH(i) = Msim(i).modelData.measures.TS.RH;
    
    MS_HD_shuff(i) = mean(Msim(i).ts_hd, 'omitnan');
    MS_RH_shuff(i) = mean(Msim(i).ts_rh, 'omitnan');
    
    VE_place_shuff(i) = mean(Msim(i).ve_place, 'omitnan');
    VE_RH_shuff(i) = mean(Msim(i).ve_rh, 'omitnan');
    
    VE_place(i) = Msim(i).modelData.measures.VE.place;
    VE_RH(i) = Msim(i).modelData.measures.VE.RH;
    
end

plotName = 'RH VE';
X = VE_RH;
Xshuf = VE_RH_shuff;
stdev = std(VE_RH, 'omitnan');
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

%% PLOT ALL CELLS
for i = 1:100
    Pnow = Msim(i).param.P;
    STnow = Msim(i).ST;
    Znow = get_MD(Pnow);
    pathPlot_hd(Pnow, STnow, Znow);
    pause; close all;
end






