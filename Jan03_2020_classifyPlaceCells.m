% January 03, 2020
% Classify place cells

count = 1;
for nn = 1:length(JZ.neurons)
    for uu = 1:length(JZ.neurons(nn))
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(nn).members(uu).P);
        ST = JZ.neurons(nn).members(uu).ST;
        deltaT = mode(diff(P(:,1)));
        
        % define a threshold based on how long the session is
        spikethresh = round((nanmax(P(:,1))/6));
        
        if length(ST) > spikethresh
            % smooth position vectors
            sigma = 2; % width of Gaussian kernel
            P = smooth_pos(P_raw, sigma);
            
            % find spatial information
            spikes = ST(ST < P(end,1) & ST > P(1,1)); % remove bad STs
            edgesT = linspace(P(1,1), P(end,1), numel(P(:,1))+1);
            binnedSpikes = histcounts(spikes,edgesT);
            
            % calculate ratemap
            nbins = 50; % 3 cm bins
            map = analyses.map(P, ST, 'smooth', 2, 'binWidth', 150/nbins);
            occupancy = map.time; % time animal spent in each bin (s)
            countMap = map.count; % # of time the cell spikes in each spatial bin
            ratemap = map.z; % rate map (Hz)
            
            % firing rate in each spatial bin (Hz) 
            firing_rate = map.z;
            
            % compute probability density(proportion of time spent in each bin)
            pdx = occupancy ./ sum(occupancy, 'omitnan');
            
            % compute [overall] average firing rate (Hz)
            mean_rate = sum(binnedSpikes)/nanmax(P(:,1));
            
            % Eq(1) from Skaggs (1993) 
            SI = sum(ratemap.*log2(ratemap./mean_rate).*pdx, 'all', 'omitnan');
            
            % compute logTerm from MI equation
            logTerm = log2(firing_rate./mean_rate);
            
            % Correct for undefined values
            logTerm(countMap == 0) = 0;
            
            % Compute information value
            I = sum(firing_rate * (logTerm.*pdx)', 'all', 'omitnan');
            
            % Divide by firing rate to obtain information per spike (bits/spike) 
            Ispk = I;

            
        end
        
    end
end