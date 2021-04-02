%% Make histograms of the angular occupancy in each behavioral session

%% select session and load data
pathname = 'D:\Data\Project Data\Blackstad-OF\dset.mat';
load(pathname);

for sess = 1:86
%% (1) angular occupancy, collapsed in time
% position and head direction
P = dset(sess).P; Z = get_hd(P);
[mu, ul, ll] = circ_mean(deg2rad(Z));
avgAng(sess) = mod(rad2deg(mu),360);

end
figure(1); hold on;
aa = histcounts(avgAng,15);
bar(linspace(0,360,15),aa, 'k');
[mu, ul, ll] = circ_mean(deg2rad(avgAng));
mu = mod(rad2deg(mu),360);
ul = mod(rad2deg(ul),360); ll = mod(rad2deg(ll),360);
% xline(ul, '--r', 'LineWidth', 1); xline(ll, '--r', 'LineWidth', 1);
xline(mu, '--r','LineWidth', 1)
set(gca, 'FontSize', 15); 
box off; pbaspect([1 1 1]);
xlabel('average angle (deg)'); ylabel('sessions (count)');
xlim([-20 380]);

%%

figure(1);
plot(P(:,2), P(:,3), 'k');
set(gca, 'visible', 'off');
pbaspect([1 1 1])
t = P(:,1); x = P(:,2); y = P(:,3);
% Z = mod((Z+180),360)-180; % shift + mirror

% sampling frequency info
tpf = mode(diff(t)); % time per frame (s)

% bin arena (10x10)
nBins = 10;
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);
xCenter = (diff(xEdges)/2) + xEdges(1:end-1);
yCenter = (diff(yEdges)/2) + yEdges(1:end-1);

% bin centers vector
count = 1;
for cc = 1:length(xCenter)
    for rr = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(rr), yCenter(cc)];
        count = count+1;
    end
end

% angular bin centers
num_Z_bins = 10; % 12 deg/bin
Z_bins = linspace(0,360,num_Z_bins+1);
Z_edges = linspace(.1, 360, num_Z_bins+1);
[Z_count_all, ~, Z_idx_all] = histcounts(Z, Z_edges);

%% PLOT CIRCULAR HISTOGRAM
warning('off', 'all')
obj1 = CircHist(Z, 15);% Adjust appearance:
obj1.colorBar = 'k';  % change color of bars
obj1.avgAngH.LineStyle = '--'; % make average-angle line dashed
obj1.avgAngH.LineWidth = .10; % make average-angle line thinner
obj1.colorAvgAng = [.5 .5 .5]; % change average-angle line color
% remove offset between bars and plot-center
rl = rlim; % get current limits
obj1.setRLim([0, rl(2)]); % set lower limit to 0
% draw circle at r == 0.5 (where r == 1 would be the outer plot edge)
rl = rlim;
obj1.drawCirc((rl(2) - rl(1)) /2, '--b', 'LineWidth', 2)
obj1.scaleBarSide = 'right'; % draw rho-axis on the right side of the plot
obj1.polarAxs.ThetaZeroLocation = 'top'; % rotate the plot to have 0° on the right side
obj1.setThetaLabel('Direction', 'bottomleft'); % add label
obj1.setThetaLabel(''); % add label
% draw resultant vector r as arrow
delete(obj1.rH)
obj1.drawArrow(obj1.avgAng, obj1.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')
% Change theta- and rho-axis ticks
obj1.polarAxs.ThetaAxis.MinorTickValues = []; % remove dotted tick-lines
thetaticks(0:90:360); % change major ticks
% rticks(0:4:20); % change rho-axis tick-steps
obj1.drawScale; % update scale bar
warning('on', 'all')

stats = circ_stats(deg2rad(Z));
disp(['std: ', num2str(rad2deg(stats.std))]);
disp(['std0: ', num2str(rad2deg(stats.std0))]);

pause; clf;
end


%% plot linear histogram (boring)
% plot histogram 
figure(2)
bar(Z_edges(2:end), Z_count*tpf, 'k');
set(gca, 'FontSize', 15); box off;
xlabel('angle (deg)');
ylabel('occupancy (sec)');
pbaspect([1 1 1]);


%% angular occupancy, across space
count = 1;
rcount = 0;
for rr = 1:nBins
    rcount = rcount + 1;
    for cc = 1:nBins
        
        disp(['r = ', num2str(rr), ', c = ', num2str(cc), ', count = ', num2str(count)])
        
        % find frames in which animal occupied this spatial bin
        idx_here = find(rr == binY & cc == binX);
        time_in_bin = length(idx_here)*tpf; % occupancy (s)
        
        % spikes and angular variable Z (rad) in this spatial bin
        Z_here = Z(idx_here);
        
        here(rr,cc) = length(Z_here);
        
        Zhere{rr,cc} = Z_here;
        
        % compute occupany in each angular bin
        % linear histogram because input constrained between -pi and pi
        Z_edges = linspace(-180, 180, num_Z_bins+1);
        [Z_count_here, ~, Z_idx_here] = histcounts(Z_here, Z_edges);
        Z_idx_here(Z_idx_here==0) = NaN; 
        
        % calculate angular occupancy (s) for this spatial bin
        Z_occ_here{rr,cc} = Z_count_here .* tpf; 
        lenZocchere(rr,cc) = length(Z_count_here);
        
%         subplot(5,5,count)
%         histnow = polarhistogram(mod(Z_here,360), 30, 'Normalization', 'probability', ...
%             'FaceColor','k','FaceAlpha',.3);
%         histnow.Parent.ThetaZeroLocation = 'top';
%         histnow.Parent.ThetaTick = [0 90 180 270];
%         histnow.Parent.RTick = max(histnow.Parent.RTick);
%         histnow.Parent.Visible = 'off';
%         histnow.FaceColor = 'white';
        
        count = count + 1;
    end
end

occraw = flipud(Zhere); occhist = flipud(Z_occ_here);

savename = ['D:\Data\Project Data\Blackstad-OF\of-', num2str(sess), '.mat'];
save(savename, 'occraw', 'occhist');




