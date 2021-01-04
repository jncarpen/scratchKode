%% December 23, 2020
% Make egocentric bearing tuning curves and rank them by peak orientation

clear tcMatrix peak_orientation
count = 1;
for nn = 1:length(JZ.neurons)
    for uu = 1:length(JZ.neurons(nn).members)
        
        % pull information for this neuron
        P_raw = units2cm(JZ.neurons(nn).members(uu).P);
        t = P_raw(:,1); Fs = mode(diff(t));
        ST = JZ.neurons(nn).members(uu).ST;
        
        if length(ST) > 100
            % smooth position vectors
            clear P
            sigma = 2;
            P = smooth_pos(P_raw, sigma);
            HD = get_hd(P);

            % get angles associated with each spike
            spkhd = HD(knnsearch(t, ST));

            % make tuning curve (smooth and big bins)
            tc10 = analyses.turningCurve(spkhd, HD, Fs, 'smooth', 2, 'binWidth', 10);
            angular_bins = tc10(:,1);
            
            % find the peak rate of the coarse tuning curve
            [M, I] = max(tc10(:,2));
            peak_orientation(count) = angular_bins(I);
            
            % make more fine tuning curve
            tc = analyses.turningCurve(spkhd, HD, Fs, 'smooth', 1, 'binWidth', 3);
            tc_flip = tc(:,2)';
            
            % put tuning curve values into a matrix
            tcMatrix(count,:) = tc_flip;
            
            % add to the count
            count = count + 1;
            
        end
    end 
end


%% sort head direction tuning curves and plot in a heatmap
[B,sort_ind] = sort(peak_orientation, 'ascend');
[rows, cols] = size(tcMatrix);

sortedTC = zeros(rows, cols);

for ii = 1:length(peak_orientation)
    
    ind_now = sort_ind(ii);
    
    sortedTC(ii,:) = tcMatrix(ind_now,:);
    
end

% plot it
figure; set(gcf,'color','w');
imagesc(zscore(sortedTC, [], 2))
xticks([]); yticks([]);
colormap('jet')
c = colorbar; c.FontSize = 20; c.FontName = 'Calibri Light';
set(gca,'FontSize',20, 'FontName', 'Calibri Light', 'FontWeight', 'normal');


% show that firing rates are log-normally distributed
Y = log();
