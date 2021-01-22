% neurons in JZ that are classified as place cells (manually defined)
% January 18, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 field
unimodal = [72, 76, 77, 80, 85, ...
    104, 109, 112, 113, 115, 120, 124, 125, 128, 130, ...
    136, 180 181, 183, 184, 186, 188, 190, 191, 193, 195,...
    197, 215, 216, 218, 223, 224, 230, 231, 233, 234, 235, ...
    237, 270, 271, 272, 273, 287, 294, 311, 314, 333, 335, ...
    339, 404, 425, 426, 431, 432, 433, 435, 439, 440, 442, 444, ...
    445, 451, 452, 481, 487, 498, 501, 502, 503, 504, 506, 518, ...
    521, 522, 524, 535, 536, 540, 546, 548, 550, 566, 571, 576, ...
    579, 580, 581, 584, 589, 593, 592, 599, 606, 609, 620, 648, ...
    650, 652, 654, 659];


%% False positive rate for bursting cells
warning('on','all')
count = 1;
clear BurstFP
for i = 1:length(unimodal)
    % tell user which iteration we're on
    dispMe = strcat('running simulation ', {' '}, sprintf('%.f', i), ' of 100...');
    disp(dispMe{1,1});
    
    % grab session information
    now = unimodal(i);
    P_now = smooth_pos(units2cm(JZ.neurons(now).members(1).P),2);
    ST_now = JZ.neurons(now).members(1).ST;
    HD = get_hd(P_now);
    
    % make the ratemap
    map = analyses.map(P_now, ST_now, 'smooth', 2, 'binWidth', [3 3]);
    
    % simulate a bursty neuron
    [sim] = simulate_bursty_place(map, P_now);
    
    % save simualtion parameters
    BurstFP(i).P = P_now;
    BurstFP(i).ST = sim.ST;
    
    % perform the optimization (using a monte carlo method)
    total_iters = 100; 
    clear monte error_for_comparison
    for optim_iter = 1:total_iters 
        % choose some initial conditions randomly
        initial = choose_initial_conditions(P_now);
        
        % run the model 
        [model] = modelMe(P_now, sim.ST, HD, initial);

        monte(optim_iter).model = model; % save all the runs to be compared
        error_for_comparison(optim_iter) = model.err;
    end
    
    % find run that yieled smallest error value
    [~, errorValsMin] = nanmin(error_for_comparison);

    % save the run that gives the global min
    % optimSig(nn).members(uu).modelReal = monte(errorValsMin).model;
    modelData = monte(errorValsMin).model;
    
    % save model (from simulated data)
    BurstFP(i).modelData = modelData;
    
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
        initial = choose_initial_conditions(P_now);

        % run the model on the shuffled data
        [model_shuffled] = modelMe(P_now, sim.ST, shuffled_hd, initial);

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
    BurstFP(i).g = param_g;
    BurstFP(i).thetaP = param_thetaP;
    BurstFP(i).xref = param_xref;
    BurstFP(i).yref = param_yref;
    BurstFP(i).varex_place = varex_place;
    BurstFP(i).varex_rh = varex_rh;
    BurstFP(i).modstren_hd = modstren_hd;
    BurstFP(i).modstren_rh = modstren_rh;
    BurstFP(i).shuffErr = error_shuff;
    
    count = count + 1;

end
warning('on','all')


%% find percentage of cells that are significant
% @todo: implement non-parametric method

Msim = BurstFP;
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























% plot
figure(1)
pathPlot_hd(P_now, ST_now, get_hd(P_now))    
figure(2)
map = analyses.map(P_now, ST_now, 'smooth', 2, 'binWidth', 150/50); % calculate tuning curve
peakRate = nanmax(nanmax(map.z));
rate_map_title = strcat('peak fr: ', sprintf('%.2f',peakRate));
plot.colorMap(map.z)
pbaspect([1 1 1])
colormap(gca,'jet')
c2 = colorbar; c2.FontSize = 25;
set(gca,'xtick',[])
set(gca,'ytick',[])
title(rate_map_title, 'FontName', 'Calibri light', 'FontSize', 30, 'FontWeight', 'normal');
box off;
