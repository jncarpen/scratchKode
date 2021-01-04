% November 3 (cosine fit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET SESSION/UNIT INFORMATION

% session/unit of interest?
sessNum = 57; 
unitNum = 7;

% grab position data
position = pos_cm{1,sessNum};
t = position(:,1); % this needs to be in seconds
fs = mode(diff(t)); % sampling freq in seconds
x = position(:,2);
y = position(:,3);  

head_direction = hd{1,sessNum};
ST = SpikeTimes{1,sessNum}{1,unitNum};

% grab goal location for this unit
ref_point = hwCoord{1,sessNum};
xGoal = ref_point(1,1);
yGoal = ref_point(1,2);

% bin spikes
startTime = t(1); stopTime = t(end);
ST = ST(ST < stopTime & ST > startTime);
edgesT = linspace(startTime,stopTime,numel(t)+1);
binnedSpikes = histcounts(ST,edgesT);
SpkTrn = imgaussfilt(binnedSpikes, 2, 'Padding', 'replicate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET BIN CENTERS
nBins = 10; % divide the arena into 100 2D spatial bins
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);
for i = 1:length(xEdges)
    if i+1 <= length(xEdges)
        xCenter(i) = ((xEdges(i+1)-xEdges(i))/2)+xEdges(i);
    end
end

for i = 1:length(yEdges)
    if i+1 <= length(yEdges)
        yCenter(i) = ((yEdges(i+1)-yEdges(i))/2)+yEdges(i);
    end
end

count = 1;
for xx = 1:length(xCenter)
    for yy = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(xx), yCenter(yy)];
        count = count+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE VALUES
% find allocentric + egocentric 'bearing' at each timepoint (AB/EB)
alloAng = atan2d(yGoal-y, xGoal-x)+180;
egoAng = alloAng - head_direction;

% correct for negative angles (egoAng)
neg_idx = find(egoAng<0);
egoAng(neg_idx) = egoAng(neg_idx)+360;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT MODEL
% define number
numBins_HD = 10;
binWidth_deg = 360/numBins_HD;
angBins = linspace(0,360,numBins_HD);

nBins = 10;
count = 1;
for xx = 1:nBins
    for yy = 1:nBins
        for H = 1:numBins_HD
            
            R(xx, yy, H) = 1 + g(cos(theta - theta_p) - F_prime)
            
        % grab stuff in this bin
        indices = find(xx == binX & yy == binY);
        timeInBin = length(indices)*fs; % occupancy (s)

        % calculate values for current 2D spatial bin
        spikes_here = SpkTrn(indices);
        ego_here = egoAng(indices);
        
        % find egocentric bearing at time of spikes
        spikeInds = find(spikes_here > 0);
        egoSpk = ego_here(spikeInds);
        
        
        
        
        
        
        
        
    end
end


%
options = optimset('PlotFcns',@optimplotfval);
fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1.2,1];
x = fminsearch(fun,x0,options)


options = optimset('PlotFcns',@optimplotfval);
fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1.2,1];
x = fminsearch(fun,x0,options)



