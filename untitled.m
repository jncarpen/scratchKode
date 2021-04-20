
%% OCCUPANCY ACROSS SPACE
of = openfield2;
for i = 1:86
    % position info
    P = of(i).P; Z = get_hd(P);
    t = P(:,1); x = P(:,2); y = P(:,3);
    Z = mod((Z+180),360)-180; % shift + mirror

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
    num_Z_bins = 10; % 36 deg/bin
    Z_bins = linspace(0,360,num_Z_bins+1);
    Z_bin_ctrs = ((diff(Z_bins)/2) + Z_bins(1:end-1))-180;

    count = 1;
    for rr = 1:nBins
        for cc = 1:nBins
            % what bin are we in now?
            x_bin_here(rr,cc) = cc; y_bin_here(rr,cc) = rr;

            % find frames in which animal occupied this spatial bin
            idx_here = find(rr == binY & cc == binX);
            time_in_bin = length(idx_here)*tpf; % occupancy (s)

            % compute occupany in each angular bin
            Z_here = Z(idx_here);
            Z_edges = linspace(-180, 180, nBins+1);
            [Z_count_here, ~, Z_idx_here] = histcounts(Z_here, Z_edges);
            Z_idx_here(Z_idx_here==0) = NaN;
            
            stats{rr,cc} = circ_stats(deg2rad(Z_here));

            % calculate angular occupancy (s) for this spatial bin
            Z_occ_here = Z_count_here .* tpf; 
            hdcount(rr,cc,:) = Z_count_here;
            if ~isempty(Z_here)
                hskw(rr,cc) = circ_r(deg2rad(Z_here'));
            else
                hskw(rr,cc) = nan;
            end
        end
    end
        % output struct
%         ofInfo(i).skVel = std(abs(diff(P(:,2))+abs(diff(P(:,3)))));
%         ofInfo(i).hdcount = hdcount;
%         ofInfo(i).hskw = hskw;
%         ofInfo(i).skew = skewness(ofInfo(i).hdcount(:));
%         ofInfo(i).std = std(ofInfo(i).hdcount(:));
%         allStats{1,i} = stats;
        countHere{1,i} = hdcount;
        
        % quick output
%         of_skew(i) = skewness(hdcount);
%         of_std(i) = ofInfo(i).std;
%         skVel(i) = ofInfo(i).skVel;
%         duration(i) =length(P);
end

%% ALL STATS
for j = 1:86
    statsnow=allStats{1,j};
    countNow = countHere{1,j};
    count = 1;
    for r = 1:10
        for c = 1:10
            mean_now(count) = statsnow{r,c}.kurtosis;
            spatial_count = reshape(sum(countNow, 3), 100, 1);
            
            
            
            count = count + 1;
        end
    end
    skew_mean(j) = std(mean_now);
end




%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% OCCUPANCY ACROSS TIME
% warning('off', 'all')
% for i = 1:86
%     P = OF(i).P; Z = get_hd(P);
%     Z = deg2rad(mod((Z+180),360)-180);
% %     Z = Z(find(~isnan(Z)));
% %     X = P(:,2); X = X(find(~isnan(X)));
% 
% %     c = corrcoef(P(:,2),get_speed(P), 'Rows', 'Complete');
% %     hd_spd(i) = c(2);
% %     Z = Z(1:10:end);
%     % mean vector length
% %     MVL(i) = circ_r(Z);
% %     [mu(i), ~, ~] = circ_mean(Z);
% %     [thetahatVM(i), kappaVM(i)] = circ_vmpar(Z);
% %     [Svar(i), svar(i)] = circ_var(Z);
% %     pval_circsym(i) = circ_symtest(Z);
% %     [std(i), std0(i)] = circ_std(Z);
% %     [bskew(i), bskew0(i)] = circ_skewness(Z);
% %     [pval_rtest(i), z_rtest(i)] = circ_rtest(Z);
% %     [p_rao(i), U_rao(i), ~] = circ_raotest(Z);
% %     [pval_otest(i), m_otest(i)] = circ_otest(Z);
% %     [mp_moment(i),  rho_p_moment(i), mu_p_moment(i)] = circ_moment(Z);
% %     pval_medtest(i) = circ_medtest(Z);
% %     med(i) = circ_median(Z);
% %     [k_kurtosis(i), k0_kurtosis(i)] = circ_kurtosis(Z);
%     [rho_corrcl(i), pval_corrcl(i)] = circ_corrcl(Z, get_speed(P)); 
% %     slen(i) = length(P);
% end
% 
% % PLOT
% thing1 = rho_corrcl'; thing1 = [thing1(1:35); thing1(37:end)];
% thing2 = rh; thing2 = [thing2(1:35); thing2(37:end)];
% thing1 = [ones(length(thing1),1), thing1];
% b1 = thing1\thing2;
% yCalc1 = thing1*b1;
% s = scatter(thing1(:,2), thing2, [30], 'k', 'filled'); hold on;
% colormap('copper')
% s.MarkerFaceAlpha = .50;
% plot(thing1(:,2),yCalc1, 'r', 'LineWidth', 1);
% set(gca, 'FontSize', 15); box off;
% Rval = corrcoef(thing1(:,2),thing2);
% title(['rho_p_moment(HDOcc) v. FP(RH); R=', num2str(Rval(2))]);







