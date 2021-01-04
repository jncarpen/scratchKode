% november 15, 2020 (sunday)

%% get information

% session number to test
sessNum = 57; 
unitNum = 7;

% grab position data
position_raw = pos_cm{1,sessNum};
t = position_raw(:,1); % this needs to be in seconds
fs = mode(diff(t)); % sampling freq in seconds
x = position_raw(:,2); x2 = position_raw(:,4); 
y = position_raw(:,3); y2 = position_raw(:,5); 

% grab goal location for this unit
ref_point = hwCoord{1,sessNum};
xGoal = ref_point(1,1);
yGoal = ref_point(1,2);

% speed
[s, ~] = get_speed(position_raw);
s = s(:,1); % grab first column

% grab spike times & spike train
ST = SpikeTimes{1,sessNum}{1,unitNum};

% speed threshold before we make the spiketrain
% Get speed at time of spike and put into vector SpikeSpeed
SpikeSpeed = interp1 (t, s, ST); %in cm/s

% Set threshold
thr_d= 4; % this is the threshold set in jercog et al. (diff for dMan)       
thr_u= 100;

% Apply threshold 
a=find(SpikeSpeed>thr_d);
b=find(SpikeSpeed<thr_u);

% make position samples NaN
x(find(s<thr_d))=NaN; x(find(s>thr_u))=NaN;
y(find(s<thr_d))=NaN; y(find(s>thr_u))=NaN;
x2(find(s<thr_d))=NaN; x2(find(s>thr_u))=NaN;
y2(find(s<thr_d))=NaN; y2(find(s>thr_u))=NaN;
position = [t, x, y, x2, y2];

% Combined threshold 
c=intersect(a,b);

% Vector with filtered spikes - based on indexing from c
SpikeSpeed_fil=ST(c);
tSpk = SpikeSpeed_fil; % spike times

% MAKE SPIKE TRAIN (bin the spikes)- this is speed-thresholded
startTime = t(1);
stopTime = t(end);

% remove spike times that are outside the range of tracking times
tSpk = tSpk(tSpk < stopTime & tSpk > startTime);

edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate

binnedSpikes = histcounts(tSpk,edgesT);

sigma = 2; % smoothing factor
SpkTrn = imgaussfilt(binnedSpikes, sigma, 'Padding', 'replicate'); % smooth spiketrain

% get head direction values
head_direction = get_hd(position);

% make sure that head direction ranges from 0-360 deg
if nanmax(head_direction) < 350
    head_direction = rad2deg(head_direction);
else
    head_direction = head_direction;
end



%% get bin center locations

% divide the arena into 100 2D spatial bins
nBins = 10;
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);
yEdges = fliplr(yEdges); % flip y-vector

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



















