function [tc] = tc_YPR(self, spikes)
%TC_YPR Tuning Curve (Yaw-Pitch-Roll)
%   INPUTS
%   'self'              opti-track struct
%   'ST'                spiketimes (s) 
%   OUTPUTS
%   'tc.yaw'            yaw tuning curve
%   'tc.pitch'          pitch tuning curve
%   'tc.roll'           roll tuning curve
%   J. Carpenter, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pull out information
t = self.timestamp; sampleRate = mode(diff(t));
x = self.position.x; y = self.position.y; 

yaw = self.direction.azimuth; 
pitch = self.direction.pitch;
roll = self.direction.roll;

for m = 1:3
    switch m
        case 1
            YPR = yaw;
        case 2
            YPR = pitch;
        case 3
            YPR = roll;
    end
    
    % find time indices when cell spikes
    idx = knnsearch(t, spikes);
    spkX = x(idx); spkY = y(idx);
    spkAng = YPR(idx);
    
    % bin everything
    nBins = 40; % 9 degree bins
    angEdges = linspace(-pi,pi,nBins+1); % rad
    [spkMap, mapAxis] = histcounts(spkAng,angEdges);
    
    % make the allMap
    for bin = 1:length(angEdges)-1
        in_this_bin = find(YPR>angEdges(bin) & YPR<=angEdges(bin+1));
        allMap(1,bin) = length(in_this_bin); 
        % allMap(1,bin) = sum(P_logical(in_this_bin)); % for speed thresh
    end

    for i = 1:length(mapAxis)
        if i+1 <= length(mapAxis)
            binCtrs(i) = ((mapAxis(i+1)-mapAxis(i))/2)+mapAxis(i);
        end
    end
    tcVals = spkMap./(allMap*sampleRate + eps); 

    % smooth tuning curve values
    tcVals = imgaussfilt(tcVals, 2, 'Padding', 'circular');
    
    % save the tuning curves in appropriate spots
    switch m
        case 1
            tc.yaw = tcVals;
        case 2
            tc.pitch = tcVals;
        case 3
            tc.roll = tcVals;
    end
    
end
    % bin centers
    tc.binCtrs = binCtrs;
end

