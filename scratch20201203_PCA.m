%% DECEMBER 03, 2020
% run PCA on the data

clear binSpk;
for sessNum = 72; %1:length(SpikeTimes)
    if ~isempty(SpikeTimes{1,sessNum})
        for unitNum = 1:length(SpikeTimes{1,sessNum})
            ST = SpikeTimes{1,sessNum}{1,unitNum};
            if ST > 200
                spiketrain = SpikeTrain{1,sessNum}{1,unitNum};
                % save binned spikes for each cell in an NxT matrix
                binSpk(unitNum, :) = spiketrain;
            end
    end
    end
end

% normalize the spikes
normspikes = zscore(binSpk);

% run the PCA
[coeff,score,latent,tsquared,explained, mu] = pca(normspikes);

% plot
figure
scatter3(score(:,1),score(:,2),score(:,3))
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

%% try k-means clustering
X = normspikes;
[idx,C] = kmeans(X, 3)

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off
