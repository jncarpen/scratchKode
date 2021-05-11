%% results for RH model
load('D:\Data\Project Data\Blackstad-OF\RH-model\dsetFinal.mat');
for sess = 1:86
    load(['D:\Data\Project Data\Blackstad-OF\RH-model\jsmod-', ...
        num2str(sess), '.mat']);
    
    % check to see if dset and unit are same length
    disp(['session # ', num2str(sess)]);
    disp(['dset: ', num2str(length(dset(sess).unit))]);
    disp(['unit: ', num2str(length(unit))]);
    
    % filter the position data
    Pnow = dset(sess).P; 
%     led = Pnow(:,2:end);
%     Pnow(:,2:end) = Pnow(:,2:end).*100;
%     speednow = get_speed(Pnow);
%     led = Pnow(:,2:end); led(speednow>80) = nan;
%     Pnow = [Pnow(:,1), led];
%     
    % save session information
    dsetRH(sess).ofSessNum = dset(sess).ofSessNum;
    dsetRH(sess).path = dset(sess).path;
    dsetRH(sess).P = Pnow;
    dsetRH(sess).perBadSample = dset(sess).perBadSample;
    
    count = 1;
    for u = 1:length(dset(sess).unit)
        % PREPROCESSING:
        % (1) RH model must have fit at least 5 bins
        binsFitByRH = mean(unit(u).out.model.Rxyh, [3], 'omitnan');
        binsFitByRH = sum(~isnan(binsFitByRH), 'all');
        % (2) firing rate must be at least 0.2 Hz
        STnow = dset(sess).unit(u).ST;
        if ~isempty(STnow)
            [spikeTrain, spikeTrain_smooth] = binSpikes(Pnow(:,1), STnow);
            instantaneous_fr = spikeTrain./mode(diff(Pnow(:,1)));
            avgRate = mean(instantaneous_fr, 'all', 'omitnan');
            % (3) must have fired >100 spikes
            numSpk = length(STnow);

            if binsFitByRH > 5
%                 if avgRate > 0.05
                    if numSpk > 100
                        dsetRH(sess).unit(count).dManId = dset(sess).unit(u).dManId;
                        dsetRH(sess).unit(count).ST = dset(sess).unit(u).ST;
                        dsetRH(sess).unit(count).sigRH = unit(u).sig;
                        dsetRH(sess).unit(count).out = unit(u).out;
                        count = count + 1;
                    end
%                 end
            end
        end
    end
end

% save('D:\Data\Project Data\Blackstad-OF\RH-model\dsetRHFinal.mat', 'dsetRH', '-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % manually check the units
    count = 1;
    for sess = 1:86
    Pnow = dsetRH(sess).P; t = Pnow(:,1);
    for u = 1:length(dsetRH(sess).unit)
        ST = dsetRH(sess).unit(u).ST;
        % remove spikes outside of interval
        ST = ST(ST < t(end) & ST > t(1));
        % bin spikes every 1/2 ms
        edgesT = t(1):0.0005:t(end);
        strain = histcounts(ST, edgesT);
        % autocorrelogram of spiketimes
%         [acf,lags,bounds] = autocorr(strain,'NumLags', 0.03/0.0005);
        rhnow = dsetRH(sess).unit(u).out.measures.TS.RH;
        if ~isempty(ST) %& dsetRH(sess).unit(u).sigRH == 1 & rhnow < .4
            fig = figure('units','normalized','outerposition',[0 0 1 1]); % make fullscreen fig
            set(gcf,'color','w');
            subplot(1,2,1)
            pathPlot_hd(Pnow, ST, get_hd(Pnow));
            title(['s# ', num2str(sess), ...
                ', unit # ', num2str(u), ' ms: ', ...
                num2str(rhnow), ', SigYN:', num2str(dsetRH(sess).unit(u).sigRH)]);
            subplot(1,2,2)
            plotMe(dsetRH(sess).unit(u).out)
            pbaspect([1 1 1]); set(gca, 'visible', 'off');
            fileBody = ['unit-', num2str(count)];
            filename = strcat('D:\egoAnalysis\April26_openfieldJSRH\', fileBody, '.png');
            saveas(fig, filename);
            count = count + 1;
            close all;
            
        end
    end
    end


    %% CALCULATE SPATIAL SPARSITY
 count = 1;
 for sess = 1:86
     s = sessNow(sess);
     P = dsetRH(sess).P;
     P = [P(:,1), P(:,2:end).*100];
     hd = get_hd(P);
%      [speed, accel] = get_speed(P);
%      cc= corrcoef(P(:,2), hd);
%      phd(sess)=cc(2);
%      cc= corrcoef(speed, hd);
%      speedhd(sess) = cc(2);
%      cc= corrcoef([accel;accel(end)], hd);
%      accelhd(sess)=cc(2);
     for u = 1:length(dsetRH(sess).unit)
         sirate(count) = dsetRH(sess).unit(u).SI_rate;
         sicontent(count) = dsetRH(sess).unit(u).SI_content;
        ST = dsetRH(sess).unit(u).ST;
        % find angles at times of spikes
        spkidx = knnsearch(P(:,1), ST);
        spk_hd = hd(spkidx);
        % mean vector length of HD tuning curve
        rAng(count) = circ_r(deg2rad(spk_hd));
        tcAng = hdTuning(P,hd,ST,10);
        N_Ang = length(tcAng);
        % find angles at times of spikes
        spkidx = knnsearch(P(:,1), ST);
        spk_hd = hd(spkidx);
        % spatial ratemap
        tcSpatial = reshape(dsetRH(sess).unit(u).out.data.rxy,100,1);
        N_Spatial = length(tcSpatial);
        % angular sparsity
        sAng(count) = 1-(1/N_Ang)*(((nansum(tcAng)).^2)/(nansum(tcAng.^2)));
        sSpatial(count) = 1-(1/N_Spatial)*(((nansum(tcSpatial)).^2)/(nansum(tcSpatial.^2)));
        count = count + 1;
     end
     
%      msnow = dsetRH(s).unit(u).out.measures.TS.RH;
%      venow = dsetRH(s).unit(u).out.measures.VE.RH*100;
%      if msnow > .7
%          pathPlot_hd(P, dsetRH(s).unit(u).ST, get_hd(P));
%          title(['MS: ', num2str(msnow), ', VE: ', num2str(venow), ', S: ', num2str(s), ...
%              ', U: ', num2str(u)]);
%          pause; clf;
%      end
     
 end
 