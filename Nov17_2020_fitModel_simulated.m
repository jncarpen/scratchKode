% november 15, 2020 (sunday)

%% get information

% session to test
sessNum = 57; 
position_raw = pos_cm{1,sessNum};

% simulate [pure] egocentric bearing cell 
ref_point = hwCoord{1,sessNum};
angle_of_interest = 90;
[ST, ~, ~] = simulate_ego_cell(position_raw, ref_point, angle_of_interest);
pathPlot_quiver(P, ST, head_direction);

% grab position data
t = position_raw(:,1); % this needs to be in seconds
fs = mode(diff(t)); % sampling freq in seconds
x = position_raw(:,2); x2 = position_raw(:,4); 
y = position_raw(:,3); y2 = position_raw(:,5); 

% grab goal location for this unit
xGoal = ref_point(1,1);
yGoal = ref_point(1,2);

% speed
[s, ~] = get_speed(position_raw);
s = s(:,1); % grab first column


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

%% visualize the cell
figure(1)
set(gcf,'color','w');
[tcVals_egoAng] = egoBearing(position, tSpk, boxCtr{1,sessNum}, ref_point, 40, "True", "deg");


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

% % plot to test
% figure
% for k=1:length(binCenters)
%     plot(binCenters(k,1), binCenters(k,2), 'o')
%     hold on;
% end


%% Compute a tuning curves for each 2D spatial bin

% define spatial bins
numBins_HD = 10; % 24 degree bins
binWidth_deg = 360/numBins_HD+1;
angBins = linspace(0,360,numBins_HD+1);

% get angular bin centers
angBinCtrs = [];
for i = 1:length(angBins)
    if i+1 <= length(angBins)
        angBinCtrs(i) = ((angBins(i+1)-angBins(i))/2)+angBins(i);
    end
end
clear rateMap rateMap_HD R_ratio R_ratio_summed pdx_hd

% iterate through every 2D spatial bin 
nBins = 10;
count = 1;
didNotPass = 0;
spatial_occupancy=zeros(nBins,nBins);
RR = zeros(nBins, nBins, numBins_HD);
backwardY = nBins:-1:1; forwardX = 1:1:nBins;
for row = 1:nBins
    for col = 1:nBins
        yNow = backwardY(row);
        xNow = forwardX(col);
        indices = find(yNow == binY & xNow == binX);
        timeInBin = length(indices)*fs; % occupancy (s)
        
        xNow_(row,col) = xNow;
        yNow_(row,col) = yNow;
        
        location_now{row,col} = [binCenters(count,1), binCenters(count,2)];
        loc3d(row,col,:) = [binCenters(count,1), binCenters(count,2)];

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
                R_xyh = r_xyh./r_xy; 
                R_ratio{row,col}(H,1) = R_xyh; % save in cell array 
                RR(row,col,H) = R_xyh;
            end
            
            R_ratio_summed(row,col) = nansum(R_ratio{row,col});

        else
           disp(strcat('bin ', sprintf('%.f', count), ' did not pass criteria...'))
           didNotPass = didNotPass + 1;

        end
        count = count+1;
    end
end

didNotPass

%% fit the model (for one 2D spatial bin)
for row = 1:nBins
    for col = 1:nBins
        if ~isempty(R_ratio{row,col})
            disp(strcat('fitting model for when x=', sprintf('%.1f',location_now{row,col}(1,1)), ' and y=', sprintf('%.1f', location_now{row,col}(1,2)), '...'))

        % normalized activity curve for this spatial bin
        R = R_ratio{row,col};

        % set initial values for parameters to be fit
        p.g = .25;
        p.thetaP = 90;
        p.xref = boxCtr{1,sessNum}(1,1);
        p.yref = boxCtr{1,sessNum}(1,2);

        % set values for stable parameters
        X = location_now{row,col}(1,1);
        Y = location_now{row,col}(1,2);
        H = angBinCtrs;

        % have the model take a first guess given the parameters we've fed in
        [firstGuess_pred, firstGuess_err] = cosFit(p,X,Y,H,R);

        % fit the model
        [bestP, bestErr(row,col)] = fit('cosErr',p,{'g','thetaP','xref','yref'},X,Y,H,R);
        
        % store fit parameters
        best_g(row,col) = bestP.g;
        best_thetaP(row,col) = bestP.thetaP;
        best_xref(row,col) = bestP.xref;
        best_yref(row,col) = bestP.yref;
        
        else
            disp('insufficient data. spatial bin will be skipped...')
          
            best_g(row,col) = NaN;
            best_thetaP(row,col) = NaN;
            best_xref(row,col) = NaN;
            best_yref(row,col) = NaN;
            
        end
    end
end

%% fit model for entire arena
% set initial values for parameters to be fit
clear H
p.g = .25;
p.thetaP = 90;
p.xref = boxCtr{1,sessNum}(1,1);
p.yref = boxCtr{1,sessNum}(1,2);

% set values for stable parameters
X = loc3d(:,:,1); %binCenters(:,1); %location_now{row,col}(1,1);
Y = loc3d(:,:,2); %binCenters(:,2); %location_now{row,col}(1,2);
for h = 1:length(angBinCtrs)
    H(:,:,h) = repmat(angBinCtrs(h),10,10);
end

R = RR;

% have the model take a first guess given the parameters we've fed in
[firstGuess_pred, firstGuess_err] = cosFit(p,X,Y,H,R);

% fit the model
[bestP, bestErr(row,col)] = fit('cosErr',p,{'g','thetaP','xref','yref'},X,Y,H,R);


%% visualize
% best fit
[bestGuess_pred, bestGuess_err] = cosFit(bestP,X,Y,H,R);

% get real tuning curve of cell
% [real_tc, ~] = egoBearing(position, tSpk, boxCtr{1,sessNum}, ref_point, 10, "False", "deg")

% get real modulation curve
for h = 1:10
    real_mod_tc(h) = mean(mean(RR(:,:,h)));
end

% make tuning curve (fit by model)
for h = 1:10
    fitted_tc(h) = mean(mean(bestGuess_pred(:,:,h)));
end

% plot tuning curve
figure
set(gcf,'color','w');
hold on;
plt_real = plot(angBinCtrs, real_mod_tc,'ko','MarkerFaceColor','k');
plt_fit = plot(angBinCtrs,fitted_tc, '--r', 'LineWidth', .90)
xlabel("head direction (deg)", 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
ylabel("activity (normalized)", 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
title("RH-angle tuning response fitted by the model", 'FontName', 'Calibri light', 'FontSize', titlefntsz, 'FontWeight', 'normal')








%%
% plot results
figure
set(gcf,'color','w');
hold on;
plot(angBinCtrs,R,'ko','MarkerFaceColor','k');
plot(angBinCtrs,firstGuess_pred, 'LineWidth', .90, 'color', [0 .5 .7])
plot(angBinCtrs,bestGuess_pred, 'LineWidth', .90, '-r')
xticks([angBinCtrs]); yticks([0, round(max(R)/2), round(max(R))])
xlabel("head direction (deg)", 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
ylabel("activity (normalized)", 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')
title("R(x,y,H) when x=56.0 and y=82.5", 'FontName', 'Calibri light', 'FontSize', titlefntsz, 'FontWeight', 'normal')
box off
l = legend('data', strcat('err(guess)=', sprintf('%.3f',firstGuess_err)), strcat('err(bestfit)=', sprintf('%.1f',bestGuess_err))); l.FontSize = 12; l.FontName = 'Calibri light';
hold off;




%% Plot
% plot properties
titlefntsz = 20;

figure(1)
set(gcf,'color','w');
set(gca, 'visible', 'off')
heatmap = imagesc(rateMap)
title("r(x,y) = Firing Rate Map", 'FontName', 'Calibri light', 'FontSize', titlefntsz, 'FontWeight', 'normal')
% alpha(0.5) 
colorbar
pbaspect([1 1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% colormap('jet')
% brighten(.5) % brighten colormap

% spatial occupancy
figure(2)
set(gcf,'color','w');
set(gca, 'visible', 'off')
heatmap = imagesc(spatial_occupancy)
title("Spatial Occupancy (s)", 'FontName', 'Calibri light', 'FontSize', titlefntsz, 'FontWeight', 'normal')
% alpha(0.5) 
colorbar
pbaspect([1 1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% colormap('jet')
% brighten(.5) % brighten colormap

figure(3)
set(gcf,'color','w');
set(gca, 'visible', 'off')
dirMod = imagesc(R_ratio_summed)
% alpha(0.5) 
colorbar
pbaspect([1 1 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% colormap('jet')
title("\Sigma(R(x,y,H)) = Directional Modulation", 'FontName', 'Calibri light', 'FontSize', titlefntsz, 'FontWeight', 'normal')
% brighten(.5)

figure(4)
set(gcf,'color','w');
pathPlot_hd(position, tSpk, get_hd(position))
goal_loc_plot = plot(xGoal, yGoal, 'o'); % Plot a line and circle markers
set(goal_loc_plot, 'markerfacecolor', 'k', 'markersize', 10)
now_plot = plot(56, 82.5, 'o'); % Plot a line and circle markers
set(now_plot, 'markerfacecolor', 'blue', 'markersize', 10)



