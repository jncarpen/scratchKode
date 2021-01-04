% check if initial conditions significantly changes the output of the model

%% use V1 data
sessNum = 2;
unitNum = 6;

% grab unit and corresponding [downsampled] data
unitdata = spikes_and_timestamps(FullNeuronBehaviorDataSet, sessNum, unitNum);

% store data in necessary variables
P = unitdata.position;
ST = unitdata.spiketimes;
hd = FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,sessNum}(:,9);

% what is the object location (in cm)
ref_point = [32, 26]; % for session 2 (first dataset)
% ref_point = [42, 17.5]; % for session 3 (first dataset)

%% use simulated data
sessNum = 57;
AOI = 90; % angle of interest
P = pos_cm{1,sessNum};
ref_point = hwCoord{1,sessNum};
[ST, ~, ~] = simulate_ego_cell(P, ref_point, AOI);

%% fit model
clear M D_pred D_initial
for i = 1:1000
    % plot
%     fig = figure('units','normalized','outerposition',[0 0 1 1]);
%     set(gcf,'color','w');

    % randomly choose intial parameters
    p = choose_initial_conditions(P);

    % fit the model
    [model] = modelMe(P, ST, p);
    M(i) = model;

    % make plot title
%     figtit = strcat('S ', sprintf('%.f', sessNum), ', U ', sprintf('%.f', unitNum), ', xref_i =', {' '}, sprintf('%.2f', p.xref(:)), ', yref_i = ', {' '}, sprintf('%.2f', p.yref(:)), ', error =', {' '}, sprintf('%.2f', d), ' cm');
%     sgtitle(figtit{1,1}, 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')

    % calculate distance between data + prediction (in cm)
    d = sqrt((ref_point(1,1)-model.bestParams.xref)^2 + (ref_point(1,2)-model.bestParams.yref)^2);
    
    D_initial(i) = sqrt((ref_point(1,1)-p.xref)^2 + (ref_point(1,2)-p.yref)^2);
    D_pred(i) = d;
%     subplot(1,2,1)
%     plot_modelDynamics(P, ST, model, ref_point)
%     plottit = strcat('xref_{initial}:', sprintf('%.1f', p.xref), ', yref_{initial}:', sprintf('%.1f', p.yref))
%     title( plottit, 'FontName', 'Calibri light', 'FontSize', 10, 'FontWeight', 'normal')
%     hold off;

%     subplot(1,2,2)
%     plot_iterations(model)
% 
%     fileBody = strcat('iter', sprintf('%.f', i));
%     filename = strcat('D:\egoAnalysis\Nov19_initialConditions_refpnt\', fileBody, '.png');
%     saveas(fig, filename);

%     close all;
end

% grab error values
for i = 1:1000
xx(i) = M(i).bestParams.xref;
yy(i) = M(i).bestParams.yref;
end
xyedges = linspace(0,60,100);
[A, xEdges, yEdges, binX, binY] = histcounts2(xx,yy,xyedges,xyedges);

imagesc(A)
colorbar
% caxis([0 10])

% calculate distance between box center and predicted points
xx=xx'; yy=yy';
[boxCtrX,boxCtrY] = getBoxCenter(P);
dist_from_ctr = sqrt((xx - boxCtrX).^2+ (yy - boxCtrY).^2);
test = [xx(find(dist_from_ctr > 100)), yy(find(dist_from_ctr > 100))];
radius = 100;

% plot
figure
plot(xx',yy', 'r.')
title('model predictions with various initial conditions of reference point')
xlabel('x'); ylabel('y');
xlim([0 60])
ylim([0 60])

% look at the results
D_initial = D_initial';

figure
plot(D_initial, zscore(err), 'k.')
xlabel('distance from object (cm)')
ylabel('error val (zscore)')
ylim([0 1])


%% compare errors + model efficiency over intial conditions

for j = 1:1000

% randomly choose intial parameters
p = choose_initial_conditions(P);
x_i(j) = p.xref;
y_i(j) = p.yref;

% fit the model
[model] = modelMe(P, ST, p);
d = sqrt((ref_point(1,1)-model.bestParams.xref)^2 + (ref_point(1,2)-model.bestParams.yref)^2);

err_cm(j) = d;
err_val(j) = model.err;

end

% remove outliers (for plotting purposes)
[~, ~, bin] = histcounts(err_cm);
err_corr = err_cm(find(bin<5))./10000000;
x_i_corr = x_i(find(bin<5));
y_i_corr = y_i(find(bin<5));


% plot the results
figure
hold on;
% pathPlot_quiver(P, ST, get_hd(P));
init_loc = scatter(x_i_corr, y_i_corr, [30], err_corr, 'filled');
colormap(jet)
colorbar
init_loc.MarkerFaceAlpha = 0.7;


%% initial conditions grid 

howMany = 10;
[initial, grid] = choose_initial_conditions_grid(P, howMany)

for row = 1:howMany
    for col = 1:howMany
        p = initial{row,col};
        [model] = modelMe(P, ST, p);
        d = sqrt((ref_point(1,1)-model.bestParams.xref)^2 + (ref_point(1,2)-model.bestParams.yref)^2);
        
        err_cm(row,col) = d;
        err_val(row,col) = model.err;
end
end

errZ = zscore(err_cm, 0, 'all');
imagesc(errZ);

% find standard deviation of the error (cm)
standard_dev = std(reshape(err_cm, 100,1));


%% change initial conditions (g)

clear p gNow d model err_cm err_val

for j = 1:1000
p = choose_initial_conditions_g()
gNow(j) = p.thetaP;

[model] = modelMe(P, ST, p);
d = sqrt((ref_point(1,1)-model.bestParams.xref)^2 + (ref_point(1,2)-model.bestParams.yref)^2);

err_cm(j) = d;
err_val(j) = model.err;
end

figure
plot(gNow, err_cm, 'k.')
xlim([0 360])
xlabel('initial value of theta')
ylabel('error (cm)')






