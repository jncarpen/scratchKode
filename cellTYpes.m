%% EVERYTHING YOU NEED TO LOOK AT A CELL TYPE

pop = hd1000; fp = hd1000FP;

for un = 1:length(fp)
    unitNow(un) = fp(un).B.unit;
end

HD_sig_NP = zeros(length(pop),1).*nan; RH_sig_NP = zeros(length(pop),1).*nan; 
HD_sig = zeros(length(pop),1).*nan; RH_sig = zeros(length(pop),1).*nan; 
for i = 1:length(fp)
    pop_idx = unitNow(i);
    % grab the modulation strengths for each unit
    MS_HD_now = pop(pop_idx).out.measures.TS.HD;
    MS_RH_now = pop(pop_idx).out.measures.TS.RH;
    
    % grab the shuffled distribution for each unit
    HD_shuf_now = fp(i).B.mshd;
    RH_shuf_now = fp(i).B.msrh;
    
    % find confidence interval for the shuffled distribution
    HD_ci = [mean(HD_shuf_now) - 2*std(HD_shuf_now), mean(HD_shuf_now) + 2*std(HD_shuf_now)];
    RH_ci = [mean(RH_shuf_now) - 2*std(RH_shuf_now), mean(RH_shuf_now) + 2*std(RH_shuf_now)];
    
    HD_sig(pop_idx) = MS_HD_now < HD_ci(1) | MS_HD_now > HD_ci(2);
    RH_sig(pop_idx) = MS_RH_now < RH_ci(1) | MS_RH_now > RH_ci(2);
    
    % sort the shuffled distribution and real values
    HD_sort = sort([HD_shuf_now, MS_HD_now], 'ascend');
    RH_sort = sort([RH_shuf_now, MS_RH_now], 'ascend');
    
    % use a non-parametric method to determine significance
    if find(HD_sort == MS_HD_now) > 950 | find(HD_sort == MS_HD_now) < 50; 
        HD_sig_NP(pop_idx) = 1; 
    else
        HD_sig_NP(pop_idx) = 0;
    end
    
    if find(RH_sort == MS_RH_now) > 950 | find(RH_sort == MS_RH_now) < 50; 
        RH_sig_NP(pop_idx) = 1;
    else
        RH_sig_NP(pop_idx) = 0;
    end
end

disp('COMPUTING MODEL SENSITIVITY MEASURES...')
disp(['HD TS (parametric):', num2str(nansum(HD_sig)), '%']);
disp(['HD TS (non-parametric):', num2str(nansum(HD_sig_NP)), '%']);
disp(['RH TS (parametric):', num2str(nansum(RH_sig)), '%']);
disp(['RH TS (non-parametric):', num2str(nansum(RH_sig_NP)), '%']);

sig_idx = find(RH_sig_NP==1);
nsig_idx = find(RH_sig_NP==0);

%% PULL OUT MODEL-FIT PARAMETERS
 % pull out model fit parameters
    g = ones(length(pop),1)*nan;
    xref = ones(length(pop),1)*nan;
    yref = ones(length(pop),1)*nan;
    theta = ones(length(pop),1)*nan;
    param_theta = ones(length(pop),1)*nan;
    param_A = ones(length(pop),1)*nan;
    param_kappa = ones(length(pop),1)*nan;
    rbar = ones(length(pop),1)*nan;
    rbar_data = ones(length(pop),1)*nan;
    rbar_shuff = ones(length(pop),1)*nan;
    rbardata_shuff = ones(length(pop),1)*nan;
    for j = 1:length(unitNow)
        i = unitNow(j);
        g(i) = pop(i).out.model.fitParams.g;
        xref(i) = pop(i).out.model.fitParams.xref*15;
        yref(i) = pop(i).out.model.fitParams.yref*15;
        theta(i) = mod(pop(i).out.model.fitParams.thetaP,360);
%         ctrx(i) = pop(i).root.ctr(1);
%         ctry(i) = pop(i).root.ctr(2);
%         sigmax(i) = pop(i).root.sigma(1);
%         sigmay(i) = pop(i).root.sigma(2);
        
%         ctrx(i) = pop(i).param.ctr(1);
%         ctry(i) = pop(i).param.ctr(2);
%         sigmax(i) = pop(i).param.sigma(1);
%         sigmay(i) = pop(i).param.sigma(2);

        param_theta(i) = mod(pop(i).param.theta,360);
        param_A(i) = pop(i).param.A;
        param_kappa(i) = pop(i).param.kappa;
%         param_x(i) = pop(i).param.rp(1);
%         param_y(i) = pop(i).param.rp(2);
% %         
        rbar(i) = pop(i).out.measures.TS.RH;
        rbar_data(i) = pop(i).out.measures.TS.HD;
        
        rbar_shuff(i) = mean(fp(j).B.msrh, 'omitnan');
        rbardata_shuff(i) = mean(fp(j).B.mshd, 'omitnan');
        g_shuff(i) =  mean(fp(j).B.g);
    end
    
% get distance between refpoint & predicted (in cm)
% d = sqrt((param_x-xref').^2 + (param_y-yref').^2);
% dt = abs(param_theta-theta); 


%% SIG (EGO) AND SIG (ALLO)
sig_in = find(xref(sig_idx)<150 & xref(sig_idx)>0 & yref(sig_idx)<150 & yref(sig_idx)>0);
sig_in = sig_idx(sig_in);
sig_out = 1:1000;
ix = union(sig_in, nsig_idx);
sig_out(ix) = nan;
sig_out = sig_out(find(~isnan(sig_out)));

%% PLOT 1D HISTOGRAMS
thing = g;
thing_shuff = g_shuff;
thing1 = thing(sig_in);
thing2 = thing(nsig_idx);
thing0 = thing(sig_out);
nBins = 20;
binrng = linspace(-5,5,nBins); %min(thing):(max(thing)-min(thing))/nBins:max(thing);  

counts_shuffle = histc(thing_shuff, binrng); hold on;
b3= bar(binrng, counts_shuffle, 'FaceColor',[.5 .5 .5]);
b3.FaceAlpha =.5;

counts0 = histc(thing0, binrng);
counts1 = histc(thing1, binrng);
counts2 = histc(thing2, binrng);
counts3 = counts0 + counts1 + counts2;
counts4 = counts0 + counts1;
figure(1)
b1 = bar(binrng, counts3, 'k'); %UC
hold on
b3 = bar(binrng, counts4, 'FaceColor', 'b'); % ALLO
b2 = bar(binrng, counts1, 'r'); % EGO
legend('Data 1', 'Data 2');
set(gca, 'FontSize', 15); box off;
l = legend({'shuffle','unclassified', 'RHT-Distal', 'RHT-Local'});
% l.Location = 'northeastoutside';
xlabel('r (model)');
ylabel('# neurons');

%% CDF
histbins = 40;
thing = g; thing_shuff = g_shuff;
ll = -40; ul = 40;
h1=histogram(thing(sig_in),linspace(ll,ul,histbins),'Normalization', 'cdf'); hold on;
h2=histogram(thing(sig_out),linspace(ll,ul,histbins),'Normalization', 'cdf');
h3=histogram(thing(nsig_idx),linspace(ll,ul,histbins),'Normalization', 'cdf');
h4=histogram(thing_shuff,linspace(ll,ul,histbins),'Normalization', 'cdf');

figure(2); hold on;
binCenters = (diff(h1.BinEdges)/2) + h1.BinEdges(1:end-1);
plot(binCenters,h1.Values,'LineWidth',1);
plot(binCenters,h2.Values,'LineWidth',1);
plot(binCenters,h3.Values,'LineWidth',1);
plot(binCenters,h4.Values,'LineWidth',1);
l = legend({'ego', 'allo', 'other', 'shuffle'});
l.Location = 'northeastoutside';
set(gca, 'FontSize', 15); box off;
xlabel('r (model)'); ylabel('cdf');


%% PLOT 2D HISTOGRAMS
customColor = flipud(hot);
customColor = customColor(10:length(customColor)-10,:);
nBins = 20;
[N, xEdges, yEdges, binX, binY] = histcounts2(ctrMassX,ctrMassY,nBins);
N(N==0)=nan;
imagescwithnan(N,customColor,[1 1 1]); colorbar
set(gca, 'FontSize', 15); box off;    

for i = 1:length(pop)
    t = pop(i).param.P(:,1); % timestamps
    tpf = mode(diff(t)); % time per frame
    ST = pop(i).ST; % spiketimes
    [spikeTrain, ~] = binSpikes(t, ST); % spiketrain (counts)
    inst_fr = imgaussfilt(spikeTrain./tpf,50); 
    avg_fr(i) = mean(inst_fr);
    max_fr(i) = max(inst_fr);
end

%% REFERENCE POINT SCATTER PLOTS
xrefAdj = zeros(1,length(pop)).*nan;
yrefAdj = zeros(1,length(pop)).*nan;
ctrX = 0; ctrY = 0;
for i = 1:length(pop)
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
s3 = scatter(xrefAdj(nsig_idx), yrefAdj(nsig_idx), [20], 'k', 'filled');
s1.MarkerFaceAlpha =.75; s2.MarkerFaceAlpha =.75; s3.MarkerFaceAlpha =.50;
set(gca, 'FontSize', 15); box off;
xline(75, '--k'); xline(-75, '--k');
yline(75, '--k'); yline(-75, '--k');
xlabel('model-predicted RP (x)');
ylabel('model-predicted RP (y)');
l = legend({'ego', 'allo', 'unclassified'});
l.Location = 'northeastoutside';

%% MEAN/PEAK FIRING RATE
for i = 1:length(pop)
    t = pop(i).param.P(:,1); % timestamps
    tpf = mode(diff(t)); % time per frame
    ST = pop(i).ST; % spiketimes
    [spikeTrain, ~] = binSpikes(t, ST); % spiketrain (counts)
    inst_fr = imgaussfilt(spikeTrain./tpf,50); 
    avg_fr(i) = mean(inst_fr);
    max_fr(i) = max(inst_fr);
end

%% PLACE CELL INFORMATION
ctrMassX = zeros(length(pop),1)*nan;
ctrMassY = zeros(length(pop),1)*nan;
mean_infield_rate = zeros(length(pop),1)*nan;
totalArea = zeros(length(pop),1)*nan;
for i = 1:length(pop)
    P = pop(i).param.P;
    ST = pop(i).ST;
    map = analyses.map(P,ST);
    [fieldMap, fields] = analyses.placefield(map, 'threshold', 0.2, 'minBins', 20, 'minPeak', 0.3);
    borderCov(i) = analyses.borderCoverage(fields);
    % find the largest field
    if ~isempty(fields)
        for fieldNow = 1:length(fields)
            fieldArea(fieldNow) = fields(fieldNow).area;
        end
        totalArea(i) = sum(fieldArea);
        [~, bigField] = nanmax(fieldArea);
        ctrMassX(i) = fields(bigField).x.*2.5;
        ctrMassY(i) = fields(bigField).y.*2.5;
        mean_infield_rate(i) = fields(bigField).meanRate;
        clear fieldArea
    end
end


%% VISUALIZE SOME CELLS
for j = 1:length(sig_in)
    i = now; %sig_in(j);
    rhNow = pop(i).out.measures.TS.RH;
%     if rhNow > .65
        pathPlot_hd(pop(i).param.P, pop(i).ST, pop(i).param.Z); hold on;
        scatter(xref(i), yref(i), [40], 'k', 'filled');
        title(['Unit ', num2str(i), ' r= ', num2str(rhNow)]);
        pause; clf;
%     end
end

%% FIND FALSE POSITIVE PER SESSION
xlist = [1:86]; sesslist = repmat(xlist,1,12);
sesslist = sesslist(1:1000);

for i = 1:max(xlist)
    sessIdxNow = find(sesslist==i);
    totalRepeats(i) = length(sessIdxNow);
    numberEgo(i) = length(intersect(sig_in, sessIdxNow));
    numberAllo(i) = length(intersect(sig_out, sessIdxNow));
    numberUC(i) = length(intersect(nsig_idx, sessIdxNow));
end

propEgo = (numberEgo./totalRepeats);
propAllo = (numberAllo./totalRepeats);
propUC = (numberUC./totalRepeats);

% PLOT
counts1 = histc(propUC, linspace(0,1,15));
counts2 = histc(propEgo, linspace(0,1,15));
counts3 = counts1 + counts2;
b1=bar(linspace(0,1,15), counts3, 'm'); hold on
b2=bar(linspace(0,1,15), counts1, 'b');
set(gca, 'FontSize', 15); box off;
xlabel('misclassification rate');
ylabel('# sessions');
l = legend({'false positive (ego or allo)'});


%% False positive figure
% [P P+E P+H H E]
countsEgo = [225/1000 289/1000 119/1000 0/1000 851/1000];
countsAllo = [21/1000 218/1000 412/1000 1000/1000 6/1000];
countsAll = countsEgo + countsAllo;
b1=bar(1:5, countsAll, 'k'); hold on;
b2=bar(1:5, countsEgo, 'r');
set(gca, 'FontSize', 15); box off;
ylabel('proportion tuned to RH-angle');
xlabel('cell type');
l = legend({'RH-tuned (allo)', 'RH-tuned (ego)'});

%% Scatter w/ regression line
thing1 = param_theta; thing2 = theta;
% thing1 = thing1(~isnan(thing1)); thing2 = thing2(~isnan(thing2));
thing1 = [ones(length(thing1),1), thing1];
b1 = thing1\thing2;
yCalc1 = thing1*b1;
t1 = thing1(:,2);
s1 = scatter(t1(nsig_idx), thing2(nsig_idx), [20], 'k', 'filled'); hold on;
s2 = scatter(t1(sig_out), thing2(sig_out), [20], 'b', 'filled');
s3 = scatter(t1(sig_in), thing2(sig_in), [20], 'r', 'filled'); 
s1.MarkerFaceAlpha = .50; s2.MarkerFaceAlpha = .50; s3.MarkerFaceAlpha = .50;
plot(t1,yCalc1, 'g', 'LineWidth', 1);
set(gca, 'FontSize', 15); box off;
[Rval, pval] = corrcoef(thing1(:,2),thing2);
title(['R= ', num2str(Rval(2)), ', p = ', num2str(pval(2))]);


