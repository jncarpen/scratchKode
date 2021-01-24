% Jan 22, 2020.
% Pablo Jercog sent me some code to fix the model. This script is where I
% will fix my code and check to see where I went wrong.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETUP

% inputs
P = Cell1_Trajectory;
ST = Cell1_Spikes(:,1);
t = P(:,1); 
y = P(:,3);
x = P(:,2);

% sampling frequency info
tpf = mode(diff(t)); % time per frame (s)
fps = 1/tpf; % frames/sec (Hz)

% compute MD
clear MD
for ts = 1:length(t)-1  
    if ts == 1
        MD(1) = NaN;
    elseif ts == length(t)-1
        MD(end-1:end+2) = NaN;
    else
        % -pi to pi
        MD(ts) = atan2( y(ts+1) - y(ts-1), x(ts+1) - x(ts-1));
    end
end

% choose an angular variable (HD or MD);
ang_var = MD;

% remove spikes outside of viable range
startTime = t(1); stopTime = t(end);
ST = ST(ST < stopTime & ST > startTime);

% raw spike train
edgesT = linspace(startTime,stopTime,numel(t)+1);
SpkTrn = histcounts(ST,edgesT);

% bin arena (10x10)
nBins = 10;
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);

% get bin centers
for i = 1:length(xEdges)
    if i+1 <= length(xEdges)
        xCenter(i) = ((xEdges(i+1)-xEdges(i))/2)+xEdges(i);
        yCenter(i) = ((yEdges(i+1)-yEdges(i))/2)+yEdges(i);
    end
end

% bin centers vector
count = 1;
for cc = 1:length(xCenter)
    for rr = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(rr), yCenter(cc)];
        count = count+1;
    end
end

% bin angular data
num_ang_bins = 10; % (pi/5 bins)
angBinWidth = 2*pi/num_ang_bins+1;
angBins = linspace(-pi,pi,num_ang_bins+1);

% get angular bin centers
angBinCtrs = [];
for i = 1:length(angBins)
    if i+1 <= length(angBins)
        angBinCtrs(i) = ((angBins(i+1)-angBins(i))/2)+angBins(i);
    end
end


%% GENERATE RATEMAPS

clear binNowX binNowY loc3d r_xy pdx_hd r_xyh R_xyh
count = 1;
for rr = 1:nBins
    for cc = 1:nBins
        
        % what bin are we in now?
        x_bin_here(rr,cc) = cc; y_bin_here(rr,cc) = rr;
        
        % coordinates (in cm) of the location of this particular bin
        loc3d(rr,cc,:) = [binCenters(count,1), binCenters(count,2)];
        
        % find frames in which animal occupied this spatial bin
        idx_here = find(rr == binY & cc == binX);
        time_in_bin = length(idx_here)*tpf; % occupancy (s)
        
        % spikes and angular variable (rad) in this spatial bin
        spikes_here = SpkTrn(idx_here); 
        ang_var_here = ang_var(idx_here); 
              
        % average rate for this spatial bin (Hz)
        r_xy_here = sum(spikes_here)/length(idx_here)*fps;
        
        % make spatial ratemap
        % @criteria: animal must have occupied each 2D spatial bin 
        % for a total of >= 1000 ms
        if time_in_bin >= 1
            r_xy(rr,cc) = r_xy_here; 
        end
        
        % compute angular occupany in each angular bin
        ang_edges = linspace(-pi, pi, nBins+1);
        [ang_count_here, ~, ang_idx_here] = histcounts(ang_var_here, angBins);
        ang_idx_here(ang_idx_here==0) = NaN; 
        
        % calculate angular occupancy (s) for this spatial bin
        ang_occ_here = ang_count_here .* tpf; 
        
        for H = 1:length(angBinCtrs)
            % amount of time animal spent in this HD bin (s)
            time_H = ang_occ_here(H);
            count_H = ang_count_here(H);
            
            % find indices when animal occupied this HD bin
            idx_H = find(ang_idx_here == H);
                
            % @criteria: animal must have occupied each angular bin for >=100 ms
            if time_H >= 0.1 
                % spiketimes in bin(x,y,H)
                spk_H = spikes_here(idx_H);

                % conditional rate, r(x,y,H)
                r_xyh_now = sum(spk_H)/(count_H*fps);
                r_xyh(rr,cc,H) = r_xyh_now;

                % conditional (normalized) rate, R(x,y,H)
                if r_xy_here == 0 && r_xyh_now == 0
                    % since occupancy criteria is met by this step,
                    % rate map should be 0 (not NaN- which you would
                    % get if r_xy_now = r_xyh_now = 0 and 0/0 = NaN).
                    R_xyh(rr,cc,H) = 0;
                else
                    % normalize by average rate in spatial bin
                    R_xyh(rr,cc,H) = r_xyh_now./r_xy_here;
                end
            end
        end  
    end
    count = count + 1; 
end


%% SEE IF THINGS LOOK RIGHT (r_xy looks wierd)
[fig, map] = plot_ratemap(P, ST);
figure; pathPlot_hd(P, ST, MD);
figure; imagesc(flipud(r_xy)); colorbar;


%% OPTIMIZATION
% randomly choose some initial conditions
p = choose_initial_conditions(P);

% min firing rate for a bin to be considered
rCutOff = .5;

%initial conditions of the 4 parameters to fit
pFitLM = [p.g; p.thetaP; p.xref; p.yref]; 

% options for fminsearch
options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8);

% perform the optimization
[pFit, fF]=fminsearch(@(pFit)aFitLMNew(pFit,R_xyh,r_xy,rCutOff,nBins),...
    pFitLM,options);

% get model-predicted firing rates for best-fit parameters
R_xyh_model = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins);


%% MODULATION STRENGTH

























