% FEB 3, 2021
% work out the bug with the variance explained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in some of pj's data
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example1_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example2_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example3_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example4_RawData.mat')
load('D:\Data\Dataset\sample data\fromPablo\Archive\Example5_RawData.mat')

% arrange dataset
Jercog(1).P = Cell1_Trajectory;
Jercog(2).P = Cell2_Trajectory;
Jercog(3).P = Cell3_Trajectory;
Jercog(4).P = Cell4_Trajectory;
Jercog(5).P = Cell5_Trajectory;

Jercog(1).ST = Cell1_Spikes;
Jercog(2).ST = Cell2_Spikes;
Jercog(3).ST = Cell3_Spikes;
Jercog(4).ST = Cell4_Spikes;
Jercog(5).ST = Cell5_Spikes;

Jercog(1).MD = get_MD(Cell1_Trajectory);
Jercog(2).MD = get_MD(Cell2_Trajectory);
Jercog(3).MD = get_MD(Cell3_Trajectory);
Jercog(4).MD = get_MD(Cell4_Trajectory);
Jercog(5).MD = get_MD(Cell5_Trajectory);

close all;
for unit = 1:5
    P = Jercog(unit).P;
    ST = Jercog(unit).ST(:,1);
    Z = rad2deg(Jercog(unit).MD);
    out = modelMe(P, ST, Z);
    thetanow = num2str(out.model.fitParams.thetaP);
    Jercog(unit).out = out;
    figure(unit)
    subplot(1,2,1)
    pathPlot_hd(P,ST,Z);
    title(thetanow);
    subplot(1,2,2)
    plotMe(out);
end













% grab some position data
whichSession = 1;
P = Jercog(whichSession).P;
Z = rad2deg(get_MD(P))';

% simulate egocentric bearing cell
param.theta = 0;
param.Z = Z;
param.P = P;
param.kappa = 5;
param.rp = [0, 0];
param.A = 10;
[sim] = simulate_ego(param);

% spiketimes
ST = sim.ST;
figure(1)
pathPlot_hd(P, ST, Z);

% optimization
out = modelMe(P, ST, Z);
figure(2)
plotMe(out);

%% histogram of jan sigurds VE values

clear VEM VEP
count = 1;
for n = 1:length(CosineFit.neuron)
    for u = 1:numel(CosineFit.neuron(n))
        VEM(count) = CosineFit.neuron(n).unit(u).out.measures.VE.RH;
        VEP(count) = CosineFit.neuron(n).unit(u).out.measures.VE.place;
        n_now(count) = n;
        u_now(count) = u;
        count = count+1;
    end
end

% 
neg = find(VEM<0);
for k = 1:length(neg)
    idxnow = neg(k);
    
end









% plot
figure(1); hold on;
histogram(ve_pablo);
histogram(VEM);
legend({'pj', 'js'});




















