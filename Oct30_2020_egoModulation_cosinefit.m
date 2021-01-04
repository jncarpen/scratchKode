% use fit to a cosine instead of MVL (test on Jan's data)

% for testing (grab session/unit number)
sessNum = 2;
position = pos_cm{1,sessNum};
ref_point = hwCoord{1,sessNum};

% grab position data
t = position(:,1); % this needs to be in seconds
fs = mode(diff(t)); % sampling freq in seconds
x = position(:,2);
y = position(:,2);    

% make sure that head direction ranges from 0-360 deg
if nanmax(tracker.hd) < 350
    head_direction = rad2deg(hd{1,sessNum});
else
    head_direction = hd{1,sessNum};
end

% grab goal location for this unit
xGoal = ref_point(1,1);
yGoal = ref_point(1,2);

% speed
[s, ~] = get_speed(position);

% grab spike times & spike train
ST = unit.t;

% speed threshold before we make the spiketrain
% Get speed at time of spike and put into vector SpikeSpeed
   SpikeSpeed = interp1 (t, s, ST);

   % Set threshold
   thr_d=0.04; % this is the threshold set in jercog et al.       
   thr_u=2;

   % Apply threshold 
   a=find(SpikeSpeed>thr_d);
   b=find(SpikeSpeed<thr_u);

   % Combined threshold 
   c=intersect(a,b);

   % Vector with filtered spikes - based on indexing from c
   SpikeSpeed_fil=ST(c);
   tSpk = SpikeSpeed_fil;

% MAKE SPIKE TRAIN (bin the spikes)- this is speed-thresholded
startTime = t(1);
stopTime = t(end);

% remove spike times that are outside the range of tracking times
tSpk = tSpk(tSpk < stopTime & tSpk > startTime);

edgesT = linspace(startTime,stopTime,numel(t)+1); % binsize is close to video frame rate

binnedSpikes = histcounts(tSpk,edgesT);

sigma = 2; % smoothing factor
SpkTrn = imgaussfilt(binnedSpikes, sigma, 'Padding', 'replicate'); % smooth spiketrain

%% get bin center locations
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

%% Calculate values
% find allocentric + egocentric 'bearing' at each timepoint (AB/EB)
alloAng = atan2d(yGoal-y, xGoal-x)+180;
egoAng = alloAng - head_direction;

% correct for negative angles (egoAng)
neg_idx = find(egoAng<0);
egoAng(neg_idx) = egoAng(neg_idx)+360;


%% Compute a tuning curves for each 2D spatial bin
numBins_HD = 40; % 9 degree bins
binWidth_deg = 360/numBins_HD;
angBins = linspace(0,360,numBins_HD);

% set up all the cell arrays youre gonna need
HD_TC = cell(10,10); EGO_TC = cell(10,10); ALLO_TC = cell(10,10);
HD_ST = cell(10,10); EGO_ST = cell(10,10); ALLO_ST = cell(10,10);
allo_angleAtPeak = []; ego_angleAtPeak = [];

% iterate through every 2D spatial bin 
nBins = 10;
count = 1;
for xx = 1:nBins
    for yy = 1:nBins
        indices = find(xx == binX & yy == binY);
        timeInBin = length(indices)*fs; % occupancy (s)

        % calculate values for current 2D spatial bin
        spikes_here = SpkTrn(indices); hd_here = head_direction(indices); 
        allo_here = alloAng(indices); ego_here = egoAng(indices);

        % find HD, AB, EB at time of spikes 
        spikeInds = find(spikes_here > 0); % since its smoothed
        angSpk = hd_here(spikeInds); 
        alloSpk = allo_here(spikeInds); egoSpk = ego_here(spikeInds);

        if ~isempty(spikeInds) && length(spikeInds) > 5 && timeInBin > 0.5
            % compute normal firing rate map [r(x,y)] for this 2D spatial bin
            rateMap_HD(xx,yy) = sum(SpkTrn(indices))./(timeInBin);

            % compute tuning curves
            % make sure that BNT is on the path (these are BNT tuning
            % curves)
            hd_tc = analyses.turningCurve(angSpk, hd_here, fs, 'smooth', 1, 'binWidth', binWidth_deg);
            allo_tc = analyses.turningCurve(alloSpk, allo_here, fs, 'smooth', 1, 'binWidth', binWidth_deg);
            ego_tc = analyses.turningCurve(egoSpk, ego_here, fs, 'smooth', 1, 'binWidth', binWidth_deg);

            % (allo): compute the orientations at the peak of the tuning curve
            allo_tc_vals = allo_tc(:,2); allo_tc_angles = allo_tc(:,1);
            [~, allo_idx] = nanmax(allo_tc_vals); 
            allo_angleAtPeak(count) = allo_tc_angles(allo_idx);


            % (ego): compute the orientations at the peak of the tuning curve
            ego_tc_vals = ego_tc(:,2); ego_tc_angles = ego_tc(:,1);
            [~, ego_idx] = nanmax(ego_tc_vals); 
            ego_angleAtPeak(count) = ego_tc_angles(ego_idx);

            % compute tuning curve statistics (also BNT)
            hd_tcStat = analyses.tcStatistics(hd_tc, binWidth_deg, 95);
            allo_tcStat = analyses.tcStatistics(allo_tc, binWidth_deg, 95);
            ego_tcStat = analyses.tcStatistics(ego_tc, binWidth_deg, 95);


            % put everything into cell arrays (10x10)
            HD_TC{xx,yy} = hd_tc; ALLO_TC{xx,yy} = allo_tc; EGO_TC{xx,yy} = ego_tc;
            HD_ST{xx,yy} = hd_tcStat; ALLO_ST{xx,yy} = allo_tcStat; EGO_ST{xx,yy} = ego_tcStat;

            % get values for peak direction
            hd_PD(count)=hd_tcStat.peakDirection;
            allo_PD(count)=allo_tcStat.peakDirection;
            ego_PD(count)=ego_tcStat.peakDirection;

            % get values for mean direction
            hd_MD(count)=hd_tcStat.mean;
            allo_MD(count)=allo_tcStat.mean;
            ego_MD(count)=ego_tcStat.mean;

            % get values for MVL
            hd_MVL(count)=hd_tcStat.r;
            allo_MVL(count)=allo_tcStat.r;
            ego_MVL(count)=ego_tcStat.r;

         else

            hd_PD(count)= NaN; allo_PD(count)= NaN; ego_PD(count)= NaN;
            hd_MD(count)= NaN; allo_MD(count)= NaN; ego_MD(count)= NaN;
            hd_MVL(count)= NaN; allo_MVL(count)=NaN; ego_MVL(count)= NaN;

            HD_TC{xx,yy} = NaN; ALLO_TC{xx,yy} = NaN; EGO_TC{xx,yy} = NaN;
            HD_ST{xx,yy} = NaN; ALLO_ST{xx,yy} = NaN; EGO_ST{xx,yy} = NaN;

        end
        count = count+1;
    end
end

%% make the quiverplot

% choose a theta vector
theta = ego_PD'; % peak firing rate

% scale values by MVL
scale = ego_MVL'; % mean vector length


% scale the scaling factor (make arrows bigger)
fac = .25;
scale = scale.*fac;

%%%%

% % test 90s vector
% theta = ones(length(theta),1)*90;
% scale = ones(length(theta),1)*6;

% % calculate an offset value
% offset = atan2d(yGoal-binCenters(:,2), xGoal-binCenters(:,1))+180;
% 
% % shift theta values (circularly)
% theta_shifted = theta+offset;
% neg_idx = find(theta_shifted<0);
% theta_shifted(neg_idx) = theta_shifted(neg_idx)+360;

theta_shifted = theta;

% Data is organized as (x, y, theta in degrees)
data = [];
data = [binCenters(:,1), binCenters(:,2), theta_shifted];

% find [unscaled] vector components (u,v)
% Based on equations: x = x0 + r*cos(theta), y = y0 + r*sin(theta)
u = cos(data(:,3) * pi/180); 
v = sin(data(:,3) * pi/180); 

% find the scaling factor (sf)
sf = abs(scale./(sqrt((u.^2)+(v.^2))));

% multiply components by scaling factor
uprime = u.*sf; %+ eps;
vprime = v.*sf; % + eps; 


%% PLOT
% format figure
% f = figure;
% set(gcf,'color','w');
% ax = gca; % current axes
% box(ax, "off")
title("test title", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
xlim([-.25 1.75]);
ylim([-.25 1.75])

hold on;
% plot theta 1 
h1 = quiver(data(:,1), data(:,2), uprime, vprime, 0); % 0 turns autoscaling off
custom_color = [.3 .5 1];
set(h1, 'Color', custom_color, 'AutoScale', 'off', 'LineWidth',.60)

% plot reference location
refPnt_plot = plot(xGoal, yGoal, 'o', 'MarkerSize', 6);
set(refPnt_plot, 'MarkerEdgeColor','r', 'markerfacecolor', 'r');

% set title
title("test title", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')

hold off;

end


%% Plot a second quiverplot (stacked)
%(can delete this later)
% theta_dos = offset;
% u2 = cos(theta_dos * pi/180); v2 = sin(theta_dos * pi/180);
% uprime2 = u2.*sf; vprime2 = v2.*sf 
% h2 = quiver(data(:,1), data(:,2), uprime2, vprime2, 0); % 0 turns autoscaling off
% set(h2, 'Color', 'red', 'AutoScale', 'off')