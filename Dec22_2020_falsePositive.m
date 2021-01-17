%% December 22, 2020
% Define the false positive rate of optimization procedure
% Using the JZ.mat data (egoHC dataset from dMan)
% J. Carpenter


%% load data/ add path 
% load("D:\Data\Dataset\dManCopyDataV1\JZ.mat") % loads in JZ.mat file
% addpath(genpath("C:\Users\17145\Documents\github_local"))
% addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001\bnt"))

count = 1;
for nn = 1:length(JZ.neurons)
    for uu = 1:length(JZ.neurons(nn))
    
    % pull information for this neuron
    P_raw = units2cm(JZ.neurons(nn).members(uu).P);
    ST = JZ.neurons(nn).members(uu).ST;
    
    if length(ST) > 100 % keep in mind that spiketimes are log-normally distributed
    
        % sampling frequency (50 samples/sec)
        Fs = mode(diff(P_raw(:,1))); 

        % smooth position vectors
        sigma = 2; % width of Gaussian kernel
        P = smooth_pos(P_raw, sigma);
    
        % get 'real' head direction values (deg)
        HD = get_hd(P);

        % perform the optimization (using a monte carlo method)
        total_iters = 100; 
        clear monte error_for_comparison
        for optim_iter = 1:total_iters 
            % choose some initial conditions randomly
            initial = choose_initial_conditions(P);
            % run the model 
            [model] = modelMe(P, ST, HD, initial);

            monte(optim_iter).model = model; % save all the runs to be compared
            error_for_comparison(optim_iter) = model.err;
        end

        % find run that yieled smallest error value
        [~, errorValsMin] = nanmin(error_for_comparison);

        % save the run that gives the global min
        % optimSig(nn).members(uu).modelReal = monte(errorValsMin).model;
        modelData = monte(errorValsMin).model;
        
        % save stuff
        M(count).modelData = modelData;
        M(count).nn = nn; M(count).uu = uu;


        %% shuffle the head direction values
        % define the number of shuffles we want
        total_shuffles = 100;
        min_shuffle = floor(60/Fs); % 30 seconds
        max_shuffle = floor(600/Fs); % 200 seconds

        % systematic shuffle (can be replaced with a pseudorandom shuffle)
        shuffle_vec = [linspace(-min_shuffle, -max_shuffle, total_shuffles/2), ...
            linspace(min_shuffle, max_shuffle, total_shuffles/2)]; 
        shuffle_vec =  floor(shuffle_vec); % make them all integers

        % clear stuff from before
        clear param_g param_thetaP param_xref param_yref ...
            varex_place varex_rh modstren_hd modstren_hd error_shuff

        % shuffle the head direction values 100x
        for shuff_num = 1:total_shuffles
            clear model_shuffled
            % grab the current shuffle value from vector
            shuff_value = shuffle_vec(shuff_num);

            % circularly shifts the elements in array A by K positions
            shuffled_hd = circshift(HD, shuff_value);

            % generate random values for parameter intial conditions
            initial = choose_initial_conditions(P);

            % run the model on the shuffled data
            [model_shuffled] = modelMe(P, ST, shuffled_hd, initial);

            % save data from each run
            % (1) best fit parameters
            param_g(shuff_num) = model_shuffled.bestParams.g;
            param_thetaP(shuff_num) = model_shuffled.bestParams.thetaP;
            param_xref(shuff_num) = model_shuffled.bestParams.xref;
            param_yref(shuff_num) = model_shuffled.bestParams.yref;
            % (2) variance explained (place & model)
            varex_place(shuff_num) = model_shuffled.varExplained.place;
            varex_rh(shuff_num) = model_shuffled.varExplained.model;
            % (3) modulation/tuning strength
            modstren_hd(shuff_num) = model_shuffled.modStrength.HD;
            modstren_rh(shuff_num) = model_shuffled.modStrength.RH;
            % (4) error
            error_shuff(shuff_num) = model_shuffled.err;
            
        end

        % do the values fall within 2 standard deviations (95% CI) of the shuffled?
        varex_place_ceil = mean(varex_place) + 2*std(varex_place);
        varex_place_floor = mean(varex_place) - 2*std(varex_place);

        varex_rh_ceil = mean(varex_rh) + 2*std(varex_rh);
        varex_rh_floor = mean(varex_rh) - 2*std(varex_rh);

        modstren_hd_ceil = mean(modstren_hd) + 2*std(modstren_hd);
        modstren_hd_floor = mean(modstren_hd) - 2*std(modstren_hd);

        modstren_rh_ceil = mean(modstren_rh) + 2*std(modstren_rh);
        modstren_rh_floor = mean(modstren_rh) - 2*std(modstren_rh);
        
        error_ceil = mean(error_shuff) + 2*std(error_shuff);
        error_floor = mean(error_shuff) - 2*std(error_shuff);
        

        % significant will give a 1, not significant will give a 0
        varex_place_sig(count) = modelData.varExplained.place > varex_place_ceil | modelData.varExplained.place < varex_place_floor;
        varex_rh_sig(count) = modelData.varExplained.model > varex_rh_ceil | modelData.varExplained.model < varex_rh_floor;
        modstren_hd_sig(count) = modelData.modStrength.HD > modstren_hd_ceil | modelData.modStrength.HD < modstren_hd_floor;
        modstren_rh_sig(count) = modelData.modStrength.RH > modstren_rh_ceil | modelData.modStrength.RH < modstren_rh_floor;
        error_sig(count) = modelData.err > error_ceil | modelData.err < error_floor;
        
        
        % save stuff
        M(count).varex_place_sig = varex_place_sig(count);
        M(count).varex_rh_sig = varex_rh_sig(count);
        M(count).modstren_hd_sig = modstren_hd_sig(count);
        M(count).modstren_rh_sig = modstren_rh_sig(count);
        M(count).error_sig = error_sig(count);
        M(count).varex_place = varex_place;
        M(count).varex_rh = varex_rh;
        M(count).modstren_hd = modstren_hd;
        M(count).modstren_rh = modstren_rh;
        M(count).shuffErr = error_shuff;
        
        count = count + 1;
    end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% False positive (with simulated cells)

warning('off','all')
clear Msim
noiseLevel = linspace(0, 0.025, 100);
count = 1;
for simcell = 1:100
    % tell user which iteration youre on
    dispMe = strcat('running simulation ', {' '}, sprintf('%.f', simcell), ' of 100...');
    disp(dispMe{1,1});
    
    % grab some position data
    whichSession = randi(length(JZ.neurons)); % randomly choose a session
    P_raw = units2cm(JZ.neurons(whichSession).members(1).P);
    P_smooth = smooth_pos(P_raw, 2);
    
    % size of arena (1 dimension)
    totalSz = 150;
    circPlaceFld = floor(rand(1)*11)+6;
    
    % set ratemap parameters
    root.peakrate = (rand(1)*15)+5;
    root.ctrX = rand(1)*totalSz;
    root.ctrY = rand(1)*totalSz;
    root.sigmaX = circPlaceFld;
    root.sigmaY = circPlaceFld;
    root.size = totalSz;
    root.bins = 50;
    root.P = P_smooth;
    
    % simulate a ratemap
    [map] = simulate_ratemap(root);
    
    % simulate place cell
    [sim] = simulate_place(map, P_smooth(:,1));
    
    % save simulation parameters
    Msim(simcell).simparam = root;
    
    % rename the variables
    ST = sim.ST;
    P = P_smooth;
    HD = get_hd(P);
    
    % save some stuff
    Msim(count).ST = ST;
   
    % sampling frequency (50 samples/sec)
    Fs = mode(diff(P(:,1))); 
    
    % get 'real' head direction values (deg)
    HD = get_hd(P);

    % perform the optimization (using a monte carlo method)
    total_iters = 100; 
    clear monte error_for_comparison
    for optim_iter = 1:total_iters 
        % choose some initial conditions randomly
        initial = choose_initial_conditions(P);
        % run the model 
        [model] = modelMe(P, ST, HD, initial);

        monte(optim_iter).model = model; % save all the runs to be compared
        error_for_comparison(optim_iter) = model.err;
    end

    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(error_for_comparison);

    % save the run that gives the global min
    % optimSig(nn).members(uu).modelReal = monte(errorValsMin).model;
    modelData = monte(errorValsMin).model;

    % save stuff
    Msim(count).modelData = modelData;
%     Msim(count).nn = nn; Msim(count).uu = uu;


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
        shuffled_hd = circshift(HD, shuff_value);

        % generate random values for parameter intial conditions
        initial = choose_initial_conditions(P);

        % run the model on the shuffled data
        [model_shuffled] = modelMe(P, ST, shuffled_hd, initial);

        % save data from each run
        % (1) best fit parameters
        param_g(shuff_num) = model_shuffled.bestParams.g;
        param_thetaP(shuff_num) = model_shuffled.bestParams.thetaP;
        param_xref(shuff_num) = model_shuffled.bestParams.xref;
        param_yref(shuff_num) = model_shuffled.bestParams.yref;
        % (2) variance explained (place & model)
        varex_place(shuff_num) = model_shuffled.varExplained.place;
        varex_rh(shuff_num) = model_shuffled.varExplained.model;
        % (3) modulation/tuning strength
        modstren_hd(shuff_num) = model_shuffled.modStrength.HD;
        modstren_rh(shuff_num) = model_shuffled.modStrength.RH;
        % (4) error
        error_shuff(shuff_num) = model_shuffled.err;

    end
    
    % save stuff
    Msim(count).varex_place = varex_place;
    Msim(count).varex_rh = varex_rh;
    Msim(count).modstren_hd = modstren_hd;
    Msim(count).modstren_rh = modstren_rh;
    Msim(count).shuffErr = error_shuff;

    count = count + 1;
end
warning('on','all')


%% plot how significance changes with noise levels
% plot LOBF
x = noiseLevel.*100;
y = MS_RH;
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, '--r', 'LineWidth', 1.5);
plot(x,y, '.', 'Color', 'k')
xlabel('percent noise'); ylabel('RH tuning strength') 
set(gca,'FontSize',20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
% grid on;


%% find percentage of cells that are significant
% @todo: implement non-parametric method
clear HD_sig RH_sig HD_sig_NP RH_sig_NP 
HD_sig_NP = zeros(100,1); RH_sig_NP = zeros(100,1); 
for i = 1:100
    % grab the modulation strengths for each unit
    MS_HD_now = Msim(i).modelData.modStrength.HD;
    MS_RH_now = Msim(i).modelData.modStrength.RH;
    
    % grab the shuffled distribution for each unit
    HD_shuf_now = Msim(i).modstren_hd;
    RH_shuf_now = Msim(i).modstren_rh;
    
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
    MS_HD(i) = Msim(i).modelData.modStrength.HD;
    MS_RH(i) = Msim(i).modelData.modStrength.RH;
    
    MS_HD_shuff(i) = mean(Msim(i).modstren_hd, 'omitnan');
    MS_RH_shuff(i) = mean(Msim(i).modstren_rh, 'omitnan');
    
    VE_place(i) = mean(Msim(i).varex_place, 'omitnan');
    VE_RH(i) = mean(Msim(1).varex_rh, 'omitnan');
    
end

plotName = 'HD Tuning Strength';
X = MS_HD;
Xshuf = MS_HD_shuff;
stdev = std(MS_HD_shuff, 'omitnan');
figure;
set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', [1 0 0]);
histogram(Xshuf, 15, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'k');
xline(mean(Xshuf) + stdev*2, ':k', 'LineWidth', 1.5); % std
xline(mean(Xshuf) - stdev*2, ':k', 'LineWidth', 1.5); % std
ylabel("frequency (count)"); xlabel("tuning strength");
title(plotName, 'FontName', 'Helvetica UI', 'FontSize', 20, 'FontWeight', 'normal');
set(gca,'FontSize',20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
l = legend('real data', 'shuffled data', '95% CI', 'Location', 'northeastoutside');
legend boxoff    
box off;

% blue = [0 .3 .8];
% red = [1 0 0];

%% how many of the significant neurons have distant locations

sigUnits = find(HD_sig == 1);
for i = 1:length(sigUnits)
    sigX(i) = Msim(sigUnits(i)).modelData.bestParams.xref;
    sigY(i) = Msim(sigUnits(i)).modelData.bestParams.yref;
end


%% pie chart
plotThis = RH_sig_NP;
plotTitle = 'significant RH tuning';
prop = sum(plotThis, 'omitnan')/length(plotThis);
pc = pie(prop, 1-prop); title(plotTitle, 'FontWeight', 'normal')
pc(2).FontSize = 15; pc(1).EdgeColor = 'white'; pc(1).FaceColor = 'black'; pc(1).FaceAlpha = .4;
set(gca,'FontSize',15, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');

%% check confidence interval
sum(X > mean(X) - stdev*2 & X < mean(X) + stdev*2); % should be around 95 (if sample is 100)

%% plot autocorrelation of head direction
% calculate and plot 
[xcf,lags,~,~] = crosscorr(HD,ST, 'NumLags',100,'NumSTD',2);
figure; set(gcf,'color','w');
stem(lags*Fs, xcf, '.k');
xlabel('lag (s)'); ylabel('xcf'); title('autocorrelation of hd')
set(gca,'FontSize',20, 'FontName', 'Calibri Light', 'FontWeight', 'normal');
box off;

%% correlation of two variables
% calculate the spiketrain from spiketimes
[spikeTrain, spikeTrain_smooth] = binSpikes(P(:,1), ST);

[s, ~] = get_speed(P);
s = s(:,1); % grab first column

% calculate and plot cross-correlation between head direction and
% spiketrain
[xcf,lags,bounds,h] = crosscorr(spikeTrain, s, 'NumLags',min_shuffle,'NumSTD',2);
close all
figure; set(gcf,'color','w');
stem(lags*Fs, xcf, '.k');
xlabel('lag (s)'); ylabel('xcf'); title('cross-correlation of speed and spiketrain')
set(gca,'FontSize',20, 'FontName', 'Calibri Light', 'FontWeight', 'normal');
box off;

%% find some of the units that pass the significance criteria
% find all indices when units were deemed significant
HD_sig_vec = find(HD_sig_NP == 1);
RH_sig_vec = find(RH_sig_NP == 1);

% choose one of these measures to look at
sig_measure = RH_sig_vec;
for j = 1:length(sig_measure)
    idx_now = sig_measure(j);
    P = Msim(j).P;
    ST = Msim(j).ST;
    
    % plot
    figure(j);
    pathPlot_hd(P, ST, get_hd(P));
end
















