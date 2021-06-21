%% ANALYZE CA1 DATASET 

% load('D:\Data\Project Data\Blackstad-OF\RH-model\dsetRHFinal.mat')

% GRAB INFORMATION FROM dsetRH struct
count = 1;
for sess = 1:length(dsetRH)
    for u = 1:length(dsetRH(sess).unit)
        isSig(count) = dsetRH(sess).unit(u).sigRH;
        numSpk(count) = length(dsetRH(sess).unit(u).ST);
        xref(count) = dsetRH(sess).unit(u).out.model.fitParams.xref*15;
        yref(count) = dsetRH(sess).unit(u).out.model.fitParams.yref*15;
        theta(count) = dsetRH(sess).unit(u).out.model.fitParams.thetaP;
        g(count) = dsetRH(sess).unit(u).out.model.fitParams.g;
        rbar(count) = dsetRH(sess).unit(u).out.measures.TS.RH;
        rbar_hd(count) = dsetRH(sess).unit(u).out.measures.TS.HD;
        mvl_omitnan = reshape(dsetRH(sess).unit(u).out.measures.MVL.RH,100,1);
        mvl_omitnan = mvl_omitnan(~isnan(mvl_omitnan));
        MVL_var(count) = var(mvl_omitnan);
        MVL_skew(count) = skewness(mvl_omitnan);
        MVL_std(count) = std(mvl_omitnan);
        sessnow(count) = sess;
        unitnow(count) = u;
        outnow{count} = dsetRH(sess).unit(u).out;
        P{count} = dsetRH(sess).P;
        ST{count} = dsetRH(sess).unit(u).ST;
        SI(count) = dsetSI(sess).unit(u).contSig;
        count = count + 1;
    end
end

%% SIG IN/ SIG OUT
sig_idx = find(isSig==1);
nsig_idx = find(isSig==0);
sig_in = find(xref(sig_idx)<150 & xref(sig_idx)>0 & yref(sig_idx)<150 & yref(sig_idx)>0);
sig_in = sig_idx(sig_in);
sig_out = 1:1000;
ix = union(sig_in, nsig_idx);
sig_out(ix) = nan;
sig_out = sig_out(find(~isnan(sig_out)));


%% PLOT 1D HISTOGRAMS
thing = g;
thing0 = thing(nsig_idx);
thing1 = thing(sig_out);
thing2 = thing(sig_in);
nBins = 20;
ll=-10; ul=10;
binrng = linspace(ll,ul,nBins); %min(thing):(max(thing)-min(thing))/nBins:max(thing);

counts0 = histc(thing0, binrng);
counts1 = histc(thing1, binrng);
counts2 = histc(thing2, binrng);
counts3 = counts0 + counts1 + counts2;
counts4 = counts1 + counts2;
figure(1)
b1 = bar(binrng, counts3, 'k'); % unclassified
hold on
b3 = bar(binrng, counts4, 'FaceColor', [.6 .6 .6]); % sig out
b2 = bar(binrng, counts2, 'b'); % sig in
legend('Data 1', 'Data 2');
set(gca, 'FontSize', 15); box off;
l = legend({'unclassified', 'RHT-D', 'RHT-L'});
% l.Location = 'northeastoutside';
xlabel('r (model)');
ylabel('# neurons');


%% PDF Normalized
thing = g;
nbins = 20; ll = -20; ul = 20;
minlen = min([length(sig_in) length(sig_out) length(nsig_idx)]);
figure; h1 = histogram(datasample(thing(nsig_idx),minlen,'replace',false),linspace(ll,ul,nbins), 'Normalization', 'pdf');
figure; h2 = histogram(datasample(thing(sig_out),minlen,'replace',false),linspace(ll,ul,nbins), 'Normalization', 'pdf');
figure; h3 = histogram(datasample(thing(sig_in),minlen,'replace',false),linspace(ll,ul,nbins), 'Normalization', 'pdf');
figure; hold on;
plot(h1.BinEdges(1:end-1),smooth(h1.Values),'LineWidth',1);
plot(h2.BinEdges(1:end-1),smooth(h2.Values),'LineWidth',1)
plot(h3.BinEdges(1:end-1),smooth(h3.Values),'LineWidth',1)
l = legend({'unclassified', 'RHT-D', 'RHT-L'});

%% CDF
histbins = 40;
thing = rbar; 
ll = -1; ul = 1;
h1=histogram(thing(sig_in),linspace(ll,ul,histbins),'Normalization', 'cdf'); hold on;
h2=histogram(thing(sig_out),linspace(ll,ul,histbins),'Normalization', 'cdf');
h3=histogram(thing(nsig_idx),linspace(ll,ul,histbins),'Normalization', 'cdf');

figure(2); hold on;
binCenters = (diff(h1.BinEdges)/2) + h1.BinEdges(1:end-1);
plot(binCenters,h1.Values,'LineWidth',1);
plot(binCenters,h2.Values,'LineWidth',1);
plot(binCenters,h3.Values,'LineWidth',1);
l = legend({'RHT-L', 'RHT-D', 'unclassified'});
l.Location = 'northeastoutside';
set(gca, 'FontSize', 15); box off;
xlabel('r (model)'); ylabel('cdf');

%% REFERENCE POINT SCATTER PLOTS
xrefAdj = zeros(1,length(rbar)).*nan;
yrefAdj = zeros(1,length(rbar)).*nan;
ctrX = 0; ctrY = 0;
for i = 1:length(rbar)
    x = xref(i)-75; y = yref(i)-75;
    if x>300 || x<-300 || y>300 || y<-300
        if x<-300 && y>-300 && y<300 % left
            xLine = -300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            yLine = (m*xLine)+b;
        elseif x>-300 && x<300 && y<-300 % bottom
            yLine = -300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            xLine = (yLine-b)/m;
        elseif x>300 && y>-300 && y<300 % right
            xLine = 300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            yLine = (m*xLine)+b;
        elseif x>-300 && x<300 && y >300 % top
            yLine = 300;
            m = (y-ctrY)/(x-ctrX);
            b = ctrY-(m*ctrX);
            xLine = (yLine-b)/m;
        end 
        xrefAdj(i) = xLine; yrefAdj(i) = yLine;
    else
        xrefAdj(i) = x; yrefAdj(i) = y;
    end
    
end

hold on;
s1 = scatter(xrefAdj(sig_in), yrefAdj(sig_in), [20], 'r', 'filled');
s2 = scatter(xrefAdj(sig_out), yrefAdj(sig_out), [20], 'b', 'filled');
s3 = scatter(xrefAdj(nsig_idx), yrefAdj(nsig_idx), [20], 'g', 'filled');
s1.MarkerFaceAlpha =.25; s2.MarkerFaceAlpha =.25; s2.MarkerFaceAlpha =.25;
set(gca, 'FontSize', 15); box off;
xline(75, '--k'); xline(-75, '--k')
yline(75, '--k'); yline(-75, '--k')
xlabel('model-predicted RP (x)')
ylabel('model-predicted RP (y)')
l = legend({'RHT-Local', 'RHT-Distal', 'unclassified'});
l.Location = 'northeastoutside';
xlim([-300 300]);
ylim([-300 300]);


%% VISUALIZE SOME CELLS
for j = 1:length(sig_out)
    i = sig_out(j);
    rhNow = rbar(i);
    pnow = P{i}; stnow = ST{i};
    fileBody = ['S', num2str(sessnow(i)), '_U', num2str(unitnow(i))];
    filename = strcat('D:\egoAnalysis\May10_CA1spikePlotsSigOut\', fileBody, '.png');
    pathPlot_hd(pnow, stnow, get_hd(pnow)); hold on;
    scatter(xref(i), yref(i), [40], 'k', 'filled');
    title(['Unit ', num2str(unitnow(i)), ', Sess ', num2str(sessnow(i)),...
        ', r= ', num2str(rhNow)]);
    fig = gcf;
    saveas(fig, filename);
    close all;
%     pause; clf;
end







