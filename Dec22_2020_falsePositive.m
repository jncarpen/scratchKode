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

clear Msim
count = 1;
for simcell = 1:100
    % grab some position data
    P_raw = units2cm(JZ.neurons(simcell).members(1).P);
    P_smooth = smooth_pos(P_raw, 2);
    
    clear simparam
    simparam.position = P_smooth;
    simparam.ctr_mass = [rand(1)*150, rand(1)*150];
    simparam.noise = 0;
    simparam.width = rand(1)*5 + 5;
    simparam.theta = rand(1)*360;
    [sim] = simulate_place_egoMod(simparam);
    
    % rename the variables
    ST = sim.spiketimes;
    P = P_smooth;
    HD = get_hd(P);
    
    % save some stuff
    Msim(count).ST = ST;
    Msim(count).P = P;
   
    % sampling frequency (50 samples/sec)
    Fs = mode(diff(P_raw(:,1))); 
    
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
    Msim(count).nn = nn; Msim(count).uu = uu;


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
        varex_place varex_rh modstren_hd modstren_hd error_shuff...
        varex_place_sig varex_rh_sig modstren_rh_sig  modstren_hd_sig...
        error_sig

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
    Msim(count).varex_place_sig = varex_place_sig(count);
    Msim(count).varex_rh_sig = varex_rh_sig(count);
    Msim(count).modstren_hd_sig = modstren_hd_sig(count);
    Msim(count).modstren_rh_sig = modstren_rh_sig(count);
    Msim(count).error_sig = error_sig(count);
    Msim(count).varex_place = varex_place;
    Msim(count).varex_rh = varex_rh;
    Msim(count).modstren_hd = modstren_hd;
    Msim(count).modstren_rh = modstren_rh;
    Msim(count).shuffErr = error_shuff;

    count = count + 1;
end


%% PLOT STUFF
% compare the results of the shuffled data with the real data

% name the variables
X = varex_rh;
Xd = modelData.varExplained.model;
plotName = 'variance explained (model)';
xname = 'variance explained';

% find the confidence interval
% SEM = std(X)/sqrt(length(X)); % standard error
% ts = tinv([0.025  0.975], length(X)-1); % t-score (student's t inverse cumulative distribution function)
% CI = mean(X) + ts*SEM; % confidence intervals 

% find +/- 2 standard deviations
stdev = std(X);

% plot
figure;
set(gcf,'color','w');
hold on;
histogram(X, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'FaceColor', [.8 .8 .8]);
xline(Xd, '--k', 'LineWidth', 1.5);
xline(mean(X) + stdev*2, ':r', 'LineWidth', 1.5); % std
xline(mean(X) - stdev*2, ':r', 'LineWidth', 1.5); % std
ylabel("frequency (count)"); xlabel(xname);
title(plotName);
set(gca,'FontSize',20, 'FontName', 'Sylfaen', 'FontWeight', 'normal');
% l = legend('shuffled data', 'real data', '95% CI', 'Location', 'northeastoutside');
% legend boxoff    
box off;

%% check confidence interval
sum(X > mean(X) - stdev*2 & X < mean(X) + stdev*2); % should be around 95 (if sample is 100)

%% plot autocorrelation of head direction
% calculate and plot 
[xcf,lags,bounds,h] = crosscorr(HD,ST, 'NumLags',100,'NumSTD',2);
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



%% Plot RH angle tuning (real v. shuffle) - J.F2B
% plot the RH angle tuning of the real data (entire population) against 
% the (mean) RH tuning of the shuffled data

% find the mean of the shuffled distribution for each unit
clear X_shuff X_data
for i = 1:100
    X_shuff(i) = mean(Msim(i).modstren_rh, 'omitnan');
    X_data(i) = Msim(i).modelData.modStrength.RH;
end

% find the mean of each of these distributions
mean_dist = [mean(X_shuff, 'omitnan'), mean(X_data, 'omitnan')];

% calculate confidence intervals
ci_X = [mean(X_shuff) - 2*std(X_shuff), mean(X_shuff) + 2*std(X_shuff)];
ci_Y = [mean(X_data) - 2*std(X_data), mean(X_data) + 2*std(X_data)];

figure; hold on;
s = scatter(X_shuff, X_data, [25], 'k', 'filled');
s.MarkerFaceAlpha = .6;
meanPoint = scatter(mean_dist(1), mean_dist(2), [60], 'r', 'filled');
meanPoint.MarkerFaceAlpha = .6;
plot(ci_X, ci_Y, '--r', 'LineWidth', 1.15)
xlabel('RH tuning strength (shuffled data)'); ylabel('RH tuning strength (real data)')
set(gca,'FontSize', 15, 'FontName', 'Leelawadee UI', 'FontWeight', 'normal');
box off;







%% SCRATCH

% % run the model on a simulated cell first
% simParams.position = P; simParams.ref_point = [75, 75];
% simParams.theta = 270; simParams.noise = 0;
% simParams.width = 10; simParams.ctr_mass = [75,75];
% [sim] = simulate_ego_cell(simParams);














