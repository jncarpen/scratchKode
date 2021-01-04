%% NOVEMBER 24, 2020- jcarpenter.
% run the model on a simulated place cell with and without modulation
% by reference-heading angle (egocentric bearing on a reference point)

% EGOCENTRIC ANGLE GUIDE (0 to 360) - 
% 0 or 360 deg:     behind animal
% 90 deg:           to animal's RH side
% 270 deg:          to animal's LH side
% 180 deg:          in front of animal

% EGOCENTRIC ANGLE GUIDE (-180 to +180) - @todo
% 0 or 360 deg:     behind animal
% 90 deg:           to animal's RH side
% 270 deg:          to animal's LH side
% 180 deg:          in front of animal

%% load in Jan Sigurd's data
load('D:\Data\Dataset\25398\25398_v2.mat')


%% simulate an egocentric bearing cell
ref_point = [67, 30]; theta = 90;
[ST, STrain, head_direction] = simulate_ego_cell(P, [67, 30], theta);

[model] = modelMe(P, ST)
xErr = abs(ref_point(1,1)-model.bestParams.xref);
yErr = abs(ref_point(1,2)-model.bestParams.yref);
thetaErr = abs(theta-model.bestParams.thetaP);

% calculate distance between predicted reference point
% and hypothesis reference point
d = sqrt((ref_point(1,1)-model.bestParams.xref)^2 + (ref_point(1,2)-model.bestParams.yref)^2);

% plot
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');
figtit = strcat('err(refX)=', sprintf('%.f', xErr), 'cm,', {' '}, ...
    'err(refY)=', sprintf('%.f', yErr), 'cm, ', {' '}, ...
    'err(theta)=', sprintf('%.f', thetaErr), 'deg');
sgtitle(figtit{1,1}, 'FontName', 'Calibri light', 'FontSize', 16, 'FontWeight', 'normal')

subplot(1,3,1)
plot_modelDynamics(P, ST, model, ref_point)

subplot(1,3,2)
plot_iterations(model)

subplot(1,3,3)
plot_vectorMod(model)



%% simulate a bunch of egocentric cells
param.position = P;
param.ctr_mass = [110, 112];
param.r = 10;
[sim] = simulate_place(param)

simdata(3).type = 'allocentric place';
simdata(3).spiketimes = sim.spiketimes;
simdata(3).position = P;
simdata(3).ctr_mass = param.ctr_mass;

% simulate place ego
param.position = P;
param.theta = 270;
param.r = 20;
param.ctr_mass = [55, 75];
param.ref_point = [55, 75];
[sim] = simulate_place_ego(param);

simdata(5).type = 'egocentric place'
simdata(5).spiketimes = sim.spiketimes;
simdata(5).position = P;
simdata(5).ctr_mass = param.ctr_mass;
simdata(5).ref_point = param.ref_point;
simdata(5).theta = param.theta;


%% simulate a cell and then run it through the model
[sim] = simulate_place_egoMod(param)
% [sim] = simulate_place(param)
P = param.position;
ST = sim.spiketimes;
[model] = modelMe(P, ST)

figure; set(gcf,'color','w');
subplot(1,3,1)
plot_vectorMod(model)
subplot(1,3,2)
plot_vectorMod_data(model)
subplot(1,3,3)
plot_vectorMod_model(model)

% plot the model's fit to the real data


%% simulate egocentric cells
% [sim] = simulate_ego_cell(param)
% [sim] = simulate_place(param)
% [sim] = simulate_place_egoMod(param)
[sim] = simulate_egoDist_cell(param)
[sim] = simulate_HD(param)
P = param.position;
ST = sim.spiketimes;

% ref_point = param.ref_point;
[model] = modelMe(P, ST);
% ref_point = [model.bestParams.xref, model.bestParams.yref];

egoRateMap(P, ST, ref_point)

figure 
plot_modelDynamics(P, ST, model, ref_point)

figure
plot_vectorMod(model)

figure
plot_vectorMod_data(model)

figure
plot_vectorMod_model(model)

figure
pathPlot_hd(P, ST, get_hd(P))
title("")

tc_shuffle(P, ST, ref_point)

tc_shuffle_poisson(P, ST, ref_point, 100)


%% plot trajectory of the animal with quiver
Ptest = P(100:200,:);
HDtest = hd_sim(100:200,:);
xTest = Ptest(:,2);
yTest = Ptest(:,3);

xTest(find(xTest == nan)) = 0;
yTest(find(yTest == nan)) = 0;
HDtest(find(HDtest == nan)) = 0;

u = cos(HDtest .* pi/180); % multiplying by pi/180 converts the values into radians
v = sin(HDtest .* pi/180);

for frame = 1:length(xTest)
    plot(xTest(frame), yTest(frame), 'o')
    hold on
end

for q = 1:length(u)
    quiver(xTest(q), yTest(q), u(q), v(q), 'k');
    hold on
end
title("allocentric head direction", 'FontSize', 25, 'FontName', 'Calibri Light', 'FontWeight', 'normal')
box off


alloTest = alloAng(100:200);
uAllo = cos(alloTest .* pi/180);
vAllo = sin(alloTest .* pi/180);

for z = 1:length(u)
    quiver(xTest(z), yTest(z), uAllo(z), vAllo(z), 'red');
    hold on
end
title("allocentric orientation to reference", 'FontSize', 25, 'FontName', 'Calibri Light', 'FontWeight', 'normal')
box off

egoTest = alloTest-HDtest;
uEgo =  cos(egoTest .* pi/180);
vEgo = sin(egoTest .* pi/180);

for z = 1:length(u)
    quiver(xTest(z), yTest(z), uEgo(z), vEgo(z), 'blue');
    hold on
end
title("egocentric orientation to reference", 'FontSize', 25, 'FontName', 'Calibri Light', 'FontWeight', 'normal')
box off

legend('HD', 'allocentric bearing', 'egocentric bearing')






