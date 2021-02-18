% February 17, 2021
% Generate MVL maps
% Need to have P, ST, and HD (0 to 360 deg)

%%
% unpack position
t = P(:,1);
x = P(:,2);
y = P(:,3);
tps = mode(diff(t));

% shift HD (where 0 deg is 'up')
HD = mod(HD + 90, 360);

% find spike indices
spkidx = knnsearch(t, ST);

% grid of locations (inside arena)
nBins = 50;
[~, xEdges, yEdges, ~, ~] = histcounts2(x,y,nBins);
xCenter = (diff(xEdges)/2) + xEdges(1:end-1);
yCenter = fliplr((diff(yEdges)/2) + yEdges(1:end-1));

% angular bins
wid = 20; % width (deg)
angbins = linspace(0,360,360/wid);

warning('off', 'all');
clear MVL mu
for r = 1:nBins
    for c = 1:nBins
        % egocentric bearing for each time point
        xref = xCenter(c);
        yref = yCenter(r);
        allo = mod(atan2d(yref-y, xref-x),360);
        ego = mod(allo-HD, 360);
        % egocentric bearing for each spike
        egospk = ego(spkidx);
        % tuning curve
        tc = analyses.turningCurve(egospk, ego, tps, 'binWidth', 20);
        % compute statistics
        [pref_angle, ~, ~] = circ_mean(deg2rad(tc(:,1)), tc(:,2));
        mvl = circ_r(deg2rad(tc(:,1)), tc(:,2));
        % save output
        MVL(r,c) = mvl;
        mu(r,c) = pref_angle;
    end
end

%% plot
% find location of max MVL
[maxValue, linearIndexesOfMaxes] = max(MVL(:));
[rowsOfMaxes, colsOfMaxes] = find(MVL == maxValue);

figure; set(gcf,'color','w'); hold on;
imagesc(MVL);
scatter(rowsOfMaxes(1), colsOfMaxes(1), 30, 'k', 'filled');
pbaspect([1 1 1])
set(gca, 'visible', 'off');






