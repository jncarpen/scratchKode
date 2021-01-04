% adding the cosine fit thing to the jercog plots that are already in dMan
% November 3, 2020 (election day!)

%% get stuff from dMan
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

count = 1;
for col = 1:length(xCenter)
    for row = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(row), yCenter(col)];
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
numBins_HD = 10; % 24 degree bins
binWidth_deg = 360/numBins_HD;
angBins = linspace(0,360,numBins_HD);

% get angular bin centers
angBinCtrs = [];
for i = 1:length(angBins)
    if i+1 <= length(angBins)
        angBinCtrs(i) = ((angBins(i+1)-angBins(i))/2)+angBins(i);
    end
end

clear rateMap rateMap_HD R_ratio R_ratio_summed pdx_hd spatial_occ

% iterate through every 2D spatial bin 
nBins = 10;
count = 1;
didNotPass = 0;
spatial_occupancy=zeros(nBins,nBins);
backwardY = nBins:-1:1; forwardX = 1:1:nBins
for row = 1:nBins
    for col = 1:nBins
        yNow = backwardY(row);
        xNow = forwardX(col);
        indices = find(yNow == binY & xNow == binX);
        timeInBin = length(indices)*fs; % occupancy (s)
        
        xNow_(row,col) = xNow;
        yNow_(row,col) = yNow;
        
        % save occupancy 
        spatial_occupancy(row,col) = timeInBin;

        % calculate values for current 2D spatial bin
        spikes_here = SpkTrn(indices); hd_here = head_direction(indices); 
        
        disp(strcat('bin ', sprintf('%.f', count), 'has ', sprintf('%.2f', sum(spikes_here)), ' spikes...'))
        
                    % find head direction at times of spikes
            spikeInds = find(spikes_here > 0); % since its smoothed
            angSpk = hd_here(spikeInds);

        if ~isempty(spikeInds) && length(spikeInds) > 5 && timeInBin > 0.5
        % if ~isempty(spikeInds) && length(spikeInds) > 5 && timeInBin > 0.5 && all(ang_occupancy > .05)
        % are at least 50 degrees covered?
        
            % compute angular occupany in each *HD* bin
            [ang_counts, edges, angBinIdx] = histcounts(mod(hd_here, 360), angBins);
            pdx_hd{row,col} = ang_counts./sum(ang_counts); % probability distribution
            ang_occupancy = ang_counts .* fs; % how many seconds the animal was in each angular bin
                        
            % compute normal firing rate map [r(x,y)] for this 2D spatial bin
            r_xy = sum(spikes_here./timeInBin); % firing rate in this spatial bin
            rateMap(row,col) = r_xy; % save for later
            
            for H = 1:length(angBinCtrs)
                idx_H = find(H == angBinIdx);
                time_H = length(idx_H)*fs;
                spk_H = spikes_here(idx_H);
                r_xyh = sum(spk_H./time_H); % firing rate of cell, conditioned on space + hd bin
                rateMap_HD{row,col}(H,1) = r_xyh;
                
                % ratio (this is what we compare with the model)
                R_xyh = r_xyh/r_xy; 
                R_ratio{row,col}(H,1) = R_xyh; % save in cell array 
            end
            
            R_ratio_summed(row,col) = nansum(R_ratio{row,col});

        else
           disp(strcat('bin ', sprintf('%.f', count), ' did not pass criteria...'))
           didNotPass = didNotPass + 1;

        end
        test(row,col) = count;
        count = count+1;
    end
end

didNotPass


%% PLOT THE RESULTS

figure(1)
set(gcf,'color','w');
set(gca, 'visible', 'off')
heatmap = imagesc(rateMap)
title("Firing Rate Map", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
% alpha(0.5) 
colorbar
% colormap('jet')
% brighten(.5) % brighten colormap

figure(2)
set(gcf,'color','w');
set(gca, 'visible', 'off')
heatmap = imagesc(fliplr(spatial_occupancy))
title("Spatial Occupancy", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
% alpha(0.5) 
colorbar
% colormap('jet')
% brighten(.5) % brighten colormap

figure(3)
set(gcf,'color','w');
set(gca, 'visible', 'off')
dirMod = imagesc(fliplr(R_ratio_summed));
% alpha(0.5) 
colorbar
% colormap('jet')
title("Directional Modulation", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
% brighten(.5)

figure(4)
set(gcf,'color','w');
pathPlot_hd(position, tSpk, get_hd(position))

figure(4)
set(gcf,'color','w');
scatter(position(:,2), position(:,3),[5],head_direction,'.')
colorbar
colormap('hsv')
title("Behavior Plot (colored by HD)", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')



%% QUIVER
% choose a theta vector
theta = ego_PD'; % peak firing rate

% scale values by MVL
scale = ego_MVL'; % mean vector length

% scale the scaling factor (make arrows bigger)
% fac = .25;
fac = 12;
scale = scale.*fac;

%%%%

% test 90s vector
% angle = 180;
% theta = ones(length(theta),1).*angle;
% scale = ones(length(theta),1)*6;

% calculate an offset value
% offset = atan2d(yGoal-binCenters(:,2), xGoal-binCenters(:,1))+180;

% shift theta values (circularly)
% theta_shifted = circ_add(theta, offset);
% neg_idx = find(theta_shifted<0);
% theta_shifted(neg_idx) = theta_shifted(neg_idx)+360;
% over_360 = find(theta_shifted > 360);
% theta_shifted(over_360) = theta_shifted(over_360) - ((floor(theta_shifted(over_360)./360)).*360);

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
f = figure;
set(gcf,'color','w');
ax = gca; % current axes
box(ax, "off")

% title('360 degrees', 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')
% xlim([-.25 1.75]);
% ylim([-.25 1.75])

hold on;
% plot theta 1 
h1 = quiver(data(:,1), data(:,2), uprime, vprime, 0); % 0 turns autoscaling off
custom_color = [.3 .5 1];
set(h1, 'Color', 'k', 'AutoScale', 'off', 'LineWidth',.60)

% plot reference location
refPnt_plot = plot(xGoal, yGoal, 'o', 'MarkerSize', 6);
set(refPnt_plot, 'MarkerEdgeColor','r', 'markerfacecolor', 'r');

% set title
title("Preferred bearing: 180 degrees", 'FontName', 'Calibri light', 'FontSize', 14, 'FontWeight', 'normal')

hold off;
 
end


%% Plot a second quiverplot (stacked)
%(can delete this later)
% theta_dos = offset;
% u2 = cos(theta_dos * pi/180); v2 = sin(theta_dos * pi/180);
% uprime2 = u2.*sf; vprime2 = v2.*sf 
% h2 = quiver(data(:,1), data(:,2), uprime2, vprime2, 0); % 0 turns autoscaling off
% set(h2, 'Color', 'red', 'AutoScale', 'off')




%% FIT THE MODEL
% for row = 1:nBins
%     for col = 1:nBins
%         for H = 1:length(angBinCtrs)
%             
%             
%         end
%         
%     end
% end