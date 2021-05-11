%% ANALYZE BEHAVIORAL CORRELATES OF FALSE POSITIVE RATE

%% LOAD DATA
% load('D:\Data\Project Data\Simulation-POP2\1000-units\place1000.mat');
% load('D:\Data\Project Data\Simulation-POP2\falsepositive\place1000FP.mat');
% load('D:\Data\Project Data\Blackstad-OF\openfield2.mat');
xlist = [1:86]; sesslist = repmat(xlist,1,12);

%% FIND INFORMATION FOR EACH SESSION
nBins = 20;
for sess = 1:length(openfield2)
    x = openfield2(sess).P(:,2);
    y = openfield2(sess).P(:,3);
    [spatialOcc{sess},xedges,yedges, binx, biny] = histcounts2(x,y,nBins);
    occ_std(sess) = var(reshape(spatialOcc{sess}, length(spatialOcc{sess}).^2, 1));

%     coverage(sess) =  analyses.arenaCoverage(openfield2(sess).P, 3, 3, [150 150]);
%     c = xcorr2(spatialOcc{sess});
%     imagesc(zscore(c)); colorbar; pause; clf;
    % border zones
%     clear zoneLeft zoneTop zoneRight zoneBottom zoneBorder zoneMiddle;
%     zoneLeft = find(binx>=0 & binx<=3);
%     zoneTop = find(biny>=17 & biny<=20);
%     zoneRight = find(binx>=17 & binx<=20);
%     zoneBottom = find(biny>=0 & biny<=3);
%     zoneBorder = unique([zoneLeft; zoneTop; zoneRight; zoneBottom]);
%     % middle zones
%     zoneMiddle = 1:length(binx);
%     zoneMiddle(zoneBorder)=nan;
%     zoneMiddle = zoneMiddle(~isnan(zoneMiddle))';
    
    
    % imagesc(zscore(spatialOcc{sess})); colorbar; pause; close all;
end

%%

