%% RUN JS DATA THROUGH GLM
addpath(genpath('/gpfs/home/dwt244/jnc'))
ctrX = 75; ctrY = 75;
for sess = 1:length(dsetRH)
    disp(['running session ', num2str(sess), ' of 86...']);
    P = dsetRH(sess).P;
    posx = P(:,2); posy = P(:,3);
    posx2 = P(:,4); posy2 = P(:,5);
    posx_c = (posx + posx2)./2;
    posy_c = (posy + posy2)./2;
    post = P(:,1); boxSize = nanmean([nanmax(posx_c) nanmax(posy_c)]);
    for u = 1:length(dsetRH(sess).unit)
        ST = dsetRH(sess).unit(u).ST;
        ST = ST(ST < post(end) & ST > post(1));
        spiketrain = histcounts(ST, linspace(post(1),post(end),numel(post)+1))';
        % check reference points
        xref = dsetRH(sess).unit(u).out.model.fitParams.xref.*15;
        yref = dsetRH(sess).unit(u).out.model.fitParams.yref.*15;
        
        if xref>300 | xref<-300 | yref>300 | yref<-300
            if x<-300 & y>-300 & y<300 % x=-300 (vert)
                xLine = -300;
                m = (yref-ctrY)/(xref-ctrX);
                b = ctrY-(m*ctrX);
                yLine = (m*xLine)+b;
            elseif x>-300 & x<300 & y<-300
                yLine = -300;
                m = (yref-ctrY)/(xref-ctrX);
                b = ctrY-(m*ctrX);
                xLine = (yLine-b)/m;
            elseif x>300 & y>-300 & y<300
                xLine = 300;
                m = (yref-ctrY)/(xref-ctrX);
                b = ctrY-(m*ctrX);
                yLine = (m*xLine)+b;
            elseif x>-300 & x<300 & y >300
                yLine = 300;
                m = (yref-ctrY)/(xref-ctrX);
                b = ctrY-(m*ctrX);
                xLine = (yLine-b)/m;
            end
            xref = xLine; yref = yLine;
            ref = [xref, yref];
        end
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
        
        dsetGLM(sess).unit(u).glmout = glmout;
    end
    %% save output
    if mod(sess,10) == 0 
        disp('saving...')
        savepath = '/gpfs/home/dwt244/jnc/output/GLM/dsetGLM.mat';
        save(savepath, 'dsetGLM', '-v7.3');
    end
end
disp('final save...')
savepath = '/gpfs/home/dwt244/jnc/output/GLM/dsetGLM.mat';
save(savepath, 'dsetGLM', '-v7.3');
