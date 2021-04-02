%% FALSE POSITIVE GLM
% March 28, 2021

%% load data and format variables
addpath(genpath('/gpfs/home/dwt244/jnc'))
savepath = '/gpfs/home/dwt244/jnc/output/GLM/glm-place.mat';
warning('off', 'all');

% load files
A = load('/gpfs/home/dwt244/jnc/data/Simulation-POP1/place100.mat');
OF = load('/gpfs/home/dwt244/jnc/data/Simulation-POP1/openfield.mat');
A = A.place100; OF = OF.openfield;

for iter = 1:100
    disp(['running iter ', num2str(iter), ' of 100...'])
    
    % data for current unit
    P = A(iter).param.P;
    ST = A(iter).ST;
    out = A(iter).out;

    % format variables:
    % (1)  x-position of left LED every 20 ms (t x 1):
    posx = P(:,2);
    % (2)  y-position of left LED every 20 ms (t x 1):
    posy = P(:,3);
    % (3)  x-position of right LED every 20 ms (t x 1):
    posx2 = P(:,4);
    % (4)  y-position of right LED every 20 ms (t x 1):
    posy2 = P(:,5);
    % (5)  x-position in middle of LEDs
    posx_c = (posx + posx2)./2;
    % (6)  y-position in middle of LEDs
    posy_c = (posy + posy2)./2;
    % (7) vector of time (seconds) at every 20 ms time bin
    post = P(:,1);
    % (8) spiketrain: vector of the # of spikes in each 20 ms time bin
    ST = ST(ST < post(end) & ST > post(1));
    spiketrain = histcounts(ST, linspace(post(1),post(end),numel(post)+1))';
    % (9) boxSize: length (in cm) of one side of the square box
    boxSize = nanmean([nanmax(posx_c) nanmax(posy_c)]);
    % (10) ref: reference point for bearing calculations (added by me)
    ref = [out.model.fitParams.xref.*15, out.model.fitParams.yref.*15];


    %% fit the model
    fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')
    fit_all_ln_models
    disp('done fitting all LN models...')

    %% find the simplest model that best describes the spike train
    fprintf('(3/5) Performing forward model selection\n')
    select_best_model
    
    %% Compute the firing-rate tuning curves
    fprintf('(4/5) Computing tuning curves\n')
    compute_all_tuning_curves

    %% plot the results
    fprintf('(5/5) Plotting performance and parameters\n') 
    plot_performance_and_parameters
    toc
    
    % package output from run
    place100(iter).glmout = glmout;
    
    if mod(iter,10) == 0
        disp(['saving iteration ', num2str(iter), ' ...']);
        save(savepath, 'place100', '-v7.3');
    end
    
    
end

    %% final save
    disp('FINAL SAVE glmout.mat...')
    save(savepath, 'place100', '-v7.3');
