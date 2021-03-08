%% Angular occupancy heatmap
% March 1, 2021

% load in some data
load('D:\Data\External Data\Blackstad-OF\dataof.mat');

sess = 5;
nbins = 10;
x = dataof(sess).P(:,2); 
y = dataof(sess).P(:,3); 
Z = get_hd(dataof(sess).P);
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nbins);
% angular bin centers
num_Z_bins = 10; % 36 deg/bin
Z_bins = linspace(0,360,num_Z_bins+1);
Z_bin_ctrs = ((diff(Z_bins)/2) + Z_bins(1:end-1));

% KL DIV
clear KL
for rr = 1:nbins
    for cc = 1:nbins
        idx_here = find(rr == binY & cc == binX);
        Z_here = Z(idx_here);
        Z_edges = linspace(0, 360, nbins+1);
        [Z_count_here, ~, Z_idx_here] = histcounts(Z_here, Z_edges);
        Z_idx_here(Z_idx_here==0) = NaN; 
        pdx = (Z_count_here./sum(Z_count_here))+eps;
        pdxu = repmat(1/10, 1, 10); % uniform distribution
        KL(rr,cc) = kldiv2(Z_bin_ctrs,pdx,pdxu);
    end
end
imagesc(KL)



