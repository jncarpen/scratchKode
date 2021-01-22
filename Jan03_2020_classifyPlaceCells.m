% January 03, 2020
%% Classify place cells

clear rate_sig cont_sig neuron unit RS CS N RS_prop CS_prop I_rate_data ...
    I_rate_content irate_array icont_array irate_distrib icont_distrib
count = 1; countNeuron = 1;
for nn = 1:length(JZ.neurons)
    % clear old variables if they exist (matters for later in the loop)
    if exist('rate_sig_neuron', 'var'); clear rate_sig_neuron; end
    if exist('cont_sig_neuron', 'var'); clear cont_sig_neuron; end
    
    for uu = 1:length(JZ.neurons(nn))
        unitCount = 1;
        
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(nn).members(uu).P);
        ST = JZ.neurons(nn).members(uu).ST;
                
        % define a threshold based on how long the session is
        %spikethresh = round((nanmax(P_raw(:,1))/6));
        
        if length(ST) > 100%spikethresh
            % smooth position vectors
            P = smooth_pos(P_raw, 2);
            t = P(:,1); 
            
            % compute information content and rate 
            [I_rate, I_content] = spatial_info(P, ST);
            
            % save information rate/content 
            I_rate_data(count) = I_rate;
            I_rate_content(count) = I_content;
            
            % get unit's spiketrain
            [spktrn, ~] = binSpikes(t, ST);
            
            % shuffled distribution
            total_shuffles = 100;
            I_rate_shuff = zeros(total_shuffles,1);
            I_content_shuff = zeros(total_shuffles,1);
            
            for shuff_num = 1:total_shuffles
                shuff_value = randi(43500) + 1500; % 30 seconds to 15 mins
                
                % will it be negative or positive?
                thedecider = rand(1);
                if thedecider > .5
                    shuff_value = shuff_value*-1;
                end
                
                % circularly shift the spike train
                spktrn_shuff = circshift(spktrn, shuff_value);
                % get spike times from the spike train
                ST_shuff = t(find(spktrn_shuff >= 1));
                % get information rate & content
                [I_rate_shuff(shuff_num), I_content_shuff(shuff_num)] = spatial_info(P, ST_shuff);
            end
            
            % save the shuffled distributions for each neuron
            irate_distrib(count) = datasample(I_rate_shuff, 1);
            icont_distrib(count) = datasample(I_content_shuff, 1);
            
            % save shuffled distributions for each unit
            irate_array{1,count} = I_rate_shuff;
            icont_array{1,count} = I_content_shuff;
            
            % get confidence intervals for shuffled distribution
            iRate_ceil = mean(I_rate_shuff, 'omitnan') + 2*std(I_rate_shuff);
            iRate_floor = mean(I_rate_shuff, 'omitnan') - 2*std(I_rate_shuff);
            iRate_mean = mean(I_rate_shuff, 'omitnan');
            
            iCont_ceil = mean(I_content_shuff, 'omitnan') + 2*std(I_content_shuff);
            iCont_floor = mean(I_content_shuff, 'omitnan') - 2*std(I_content_shuff);
            iCont_mean = mean(I_content_shuff, 'omitnan');
            
            % decide if its a place cell or not
            % significant (PC) will give a 1, not significant will give a 0
            rate_sig(count) = I_rate > iRate_ceil | I_rate < iRate_floor;
            cont_sig(count) = I_rate > iCont_ceil | I_rate < iCont_floor;
            
            % information about the cell
            neuron(count) = nn; unit(count) = uu;
            
            % save 
            rate_sig_neuron(unitCount) = rate_sig(count);
            cont_sig_neuron(unitCount) = cont_sig(count);
            
            unitCount = unitCount + 1; % increase unit count by 1
            
            count = count + 1; % increase the count by 1
        end
    end
    
    % check to see if any of the units passed the criteria for this neuron
    anyPassed = exist('cont_sig_neuron', 'var');
    
    if anyPassed == 1
    
        % if any of the units (instances in which this neuron shows up) is
        % significant, then mark the cell as significant
        RS(countNeuron) = any(rate_sig_neuron);
        CS(countNeuron) = any(cont_sig_neuron);
        N(countNeuron) = nn; % which neuron

        % what proportion of the time is the unit significant?
        RS_prop(countNeuron) = ((sum(rate_sig_neuron))/length(rate_sig_neuron));
        CS_prop(countNeuron) = ((sum(cont_sig_neuron))/length(cont_sig_neuron));

        countNeuron = countNeuron + 1;
        
    end
end



%% GENERATE PLOTS

%% (1) pie chart
prop = sum(cont_sig, 'omitnan')/length(cont_sig);
pc = pie(prop, 1-prop); title('information content (bits/spike)', 'FontWeight', 'normal')
pc(2).FontSize = 20; pc(1).EdgeColor = 'white'; pc(1).FaceColor = 'black'; pc(1).FaceAlpha = .4;
set(gca,'FontSize',20, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');


%% (2) plot the histogram of the real & shuffled population
plotName = 'information rate (bits/s)';
figure;
set(gcf,'color','w');
hold on;
histogram(I_rate_data, 25, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'FaceColor', 'r');
histogram(irate_distrib, 10, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'k');
xline(mean(irate_distrib) + stdev*2, ':k', 'LineWidth', 1.5); % std
xline(mean(irate_distrib) - stdev*2, ':k', 'LineWidth', 1.5); % std
ylabel("frequency (count)"); xlabel("info rate (bits/s)");
set(gca,'FontSize',15, 'FontName', 'Helvetica UI', 'FontWeight', 'normal');
title(plotName, 'FontName', 'Helvetica UI', 'FontSize', 20, 'FontWeight', 'normal');
l = legend('real data', 'shuffled data', '95% CI', 'Location', 'northeastoutside');
legend boxoff    
box off;
xlim([-0.05 3.05])




%% plot neurons classified as place cells

for i = 1:length(neuron)
    now = neuron(i);
    for unit = 1:length(JZ.neurons(now))
        P_now = smooth_pos(units2cm(JZ.neurons(now).members(unit).P),2);
        ST_now = JZ.neurons(now).members(unit).ST;
        
        % plot
        figure
        title_now = strcat('N', sprintf('%.f', now), ' U', sprintf('%.f', unit));
        pathPlot_hd(P_now, ST_now, get_hd(P_now));
        title(title_now)
        pause; close all;
    end
end





