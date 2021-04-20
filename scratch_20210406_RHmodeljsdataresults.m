%% results for RH model
for sess = 1:86
    load(['D:\Data\Project Data\Blackstad-OF\RH-model2\jsmod-', ...
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
                if avgRate > 0.05
                    if numSpk > 100
                        dsetRH(sess).unit(count).dManId = dset(sess).unit(u).dManId;
                        dsetRH(sess).unit(count).ST = dset(sess).unit(u).ST;
                        dsetRH(sess).unit(count).sigRH = unit(u).sig;
                        dsetRH(sess).unit(count).out = unit(u).out;
                        count = count + 1;
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%      subplot(1,2,1)
%         plotMe(unit(u).out)
%         subplot(1,2,2)
%         pathPlot_hd(dset(sess).P,dset(sess).unit(u).ST,get_hd(dset(sess).P));
%         pause; clf;
    %% manually check the units
%     Pnow = dset(sess).P; t = Pnow(:,1);
%     for u = 1:length(dset(sess).unit)
%         ST = dset(sess).unit(u).ST;
%         % remove spikes outside of interval
%         ST = ST(ST < t(end) & ST > t(1));
%         % bin spikes every 1/2 ms
%         edgesT = t(1):0.0005:t(end);
%         strain = histcounts(ST, edgesT);
%         % autocorrelogram of spiketimes
%         [acf,lags,bounds] = autocorr(strain,'NumLags', 0.03/0.0005);
%         if ~isempty(ST)
%             subplot(1,2,1)
%             pathPlot_hd(Pnow, ST, get_hd(Pnow));
%             title(['s# ', num2str(sess), ...
%                 ', unit # ', num2str(u)])
%             subplot(1,2,2)
%             plot(lags(2:end), acf(2:end), 'k');
%             pbaspect([1 1 1]); box off;
%             pause; clf;
%         end
%     end


%     %% significance
%     for u = 1:length(unit)
%         sig(u) = unit(u).sig;
%     end
% 
%     percentsig = 100*(sum(sig)./length(unit));
% 
%     disp('---------------------------------------');
%     disp(['session # ', num2str(sess), ' ----->']);
%     disp(['total # units: ', num2str(length(unit))]);
%     disp(['total # sig: ', num2str(sum(sig))]);
%     disp(['% sig: ', num2str(percentsig)]);
%     disp('---------------------------------------');
%     clear unit
% end

%% save results
% savepath = ['D:\Data\Project Data\Blackstad-OF\RH-model-results\jsmodResults-', ...
%     num2str(sess), '.mat'];


