%% calculate angular sparsity and compare with MS

% load population
load('D:\Data\Project Data\Simulation-POP2\placeego100.mat')
% load RH model fp
load('D:\Data\Project Data\Simulation-POP2\falsepositive\placeego100FP.mat');

% rename variables
pop = placeego100;
fp = placeego100FP;

%% SPARSITY CALCULATIONS
for i = 1:100
    % information for unit i
    P = pop(i).param.P;
    hd = get_hd(P);
    ST = pop(i).ST;
    
    numSpk(i) = length(ST);
    % angular ratemap
    tcAng = hdTuning(P,get_hd(P),ST, 10);
    N_Ang = length(tcAng);
    % find angles at times of spikes
    spkidx = knnsearch(P(:,1), ST);
    spk_hd = hd(spkidx);
    % mean vector length of HD tuning curve
    rAng(i) = circ_r(deg2rad(spk_hd));
    % spatial ratemap
    tcSpatial = reshape(pop(i).out.data.rxy,100,1);
    N_Spatial = length(tcSpatial);
    % angular sparsity
    sAng(i) = 1-(1/N_Ang)*(((nansum(tcAng)).^2)/(nansum(tcAng.^2)));
    sSpatial(i) = 1-(1/N_Spatial)*(((nansum(tcSpatial)).^2)/(nansum(tcSpatial.^2)));
    % modulation strength
    ms_rh(i) = pop(i).out.measures.TS.RH;
    ms_hd(i) = pop(i).out.measures.TS.HD;
    
    ve_p(i) = pop(i).out.measures.VE.place;
    ve_rh(i) = pop(i).out.measures.VE.RH;
    
    % sort the shuffled distribution and real values
        RH_sort = sort([ms_rh(i), fp(i).B.msrh], 'ascend');
        if find(RH_sort == ms_rh(i)) > 950
            RHsig(i) = 1; 
        else
            RHsig(i) = 0;
        end
        
        % same for the HD 
         HD_sort = sort([ms_hd(i), fp(i).B.mshd], 'ascend');
        if find(RH_sort == ms_hd(i)) > 950
            HDsig(i) = 1; 
        else
            HDsig(i) = 0;
        end
end


%% PLOT
thing1 = sAng;
thing1 = [ones(100,1), thing1'];
thing2 = rAng';
b1 = thing1\thing2;
yCalc1 = thing1*b1;
s = scatter(thing1(:,2), thing2, [30], RHsig, 'filled'); hold on;
colormap('copper')
s.MarkerFaceAlpha = .75;
plot(thing1(:,2),yCalc1, 'r', 'LineWidth', 1);
set(gca, 'FontSize', 15); box off;
% legend({'nsig', 'sig', 'regline'});
