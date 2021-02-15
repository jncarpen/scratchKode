function tc_shuffle_dMan(input, varargin)
%TC_SHUFFLE (dMan version)
%   As in Mimica et al. 2018
%   Output:
%   This function will output the plot of the egocentric bearing tuning
%   curve with the shuffled distribution.
%
% Alg:
% I shift the data n times (circularly, between +/-15-60 seconds) and make 
% n+1 tuning curves (of egocentric bearing). The +1 is the tuning curve of 
% the actual data, which I plot in red. Then I take the mean (at each point/bin)
% of the n shuffled tuning curves and plot, in grey, the mean +/- 2 standard deviations. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Get stuff from dMan

% parse inputs
S = dManDe.Settings();
inp = inputParser();
inp.addParameter('goalLoc', 37)
% inp.addParameter('addPeakRate',false)
inp.parse(varargin{:});
inp.KeepUnmatched = true;
P = inp.Results;

% find the unit of interest?
[unit, tracker] = dManDe.helpers.handleInput(input, varargin{:});

if ~(P.goalLoc > 0 && P.goalLoc < 38)
    warning('first varargin (goal location) can only be a number between 1 and 37')
    P.goalLoc = 37;
end  

% grab position data
t = tracker.t;
x = tracker.x;
y = tracker.y;
% can we get the x2 and y2 values?
position = [t x y];
hd = tracker.hd;

% we need to make sure that hd is in degrees and
% ranges from 0 to 360; for now I'll just display the min/max values
nanmax(hd)
nanmin(hd)

% grab goal location for this unit
xGoal = unit.recording.sessions.userData.wellLocations(P.goalLoc,1);
yGoal = unit.recording.sessions.userData.wellLocations(P.goalLoc,2);
ref_point = [xGoal, yGoal];

% get other information
ST = unit.t; % spike times (im assuming)
fs = tracker.fs; % sampling frequency?

%% I. Calculate shuffled distribution

% define the number of shuffles we want
total_shuffles = 100;

% generate random floating-point values between a & b
a = -60; b = 60;
r = (b-a).*rand(total_shuffles*10,1) + a;

% only keep values that are larger than 15(s) and smaller than -15(s)
r_thresh = r(or(r>15, r<-15));

% randomly sample from 'r_thresh' to get a random string of 'shift' values.
shiftVals = datasample(r_thresh, total_shuffles);

% clear variables
tcVals_shift = [];

for iter = 1:total_shuffles
    % how much to shift (in s) for this iteration
    shift = shiftVals(iter);
    
    % get shifted timestamps
    ST_shift = circShift_TimeStamps(position, ST, shift);
    
    % get values for tuning curve (egoBear)
    [tcVals_shift(iter,:), ~] = egoBearing(position, ST_shift, ref_point, ref_point, "False", "deg");
end

% take mean (along the columns)
mean_tc = nanmean(tcVals_shift, 1);

% take standard deviation (along the columns)
std_tc = std(tcVals_shift, [], 1);

% calculate the upper & lower bounds of y (+/- 2 stds)
yu = mean_tc + 2*std_tc;
yl = mean_tc - 2*std_tc;


%% II. Calculate 'real' tuning curve

[tcVals_real, binCtrs] = egoBearing(position, ST, ref_point, ref_point, "False", "deg");


%% III. Plot tuning curves
figure
set(gcf,'color','w');

% plot shuffled data
fill([binCtrs fliplr(binCtrs)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none')
hold all
plot(binCtrs, mean_tc, ':k')

% plot actual data
plot(binCtrs, tcVals_real, 'Color', 'r', 'LineWidth', 1.5)

% format the plot
title("Egocentric Bearing")
ylabel("fr (Hz)")
xlim([0 360])
xticks([0 90 180 270 360])
xlabel("angle (deg)")
box off

hold off;

end

