%% Check to see if the circular mean messes things up

clear g_M theta_M x_M y_M ...
    g_CM theta_CM x_CM yCM
thetaOptions = linspace(0, 360, 100);
for i = 1:100
    % grab some position data from a random session
    whichSession = randi(length(JZ.neurons));
    P = smooth_pos(units2cm(JZ.neurons(whichSession).members(1).P),3);
    
%     % set parameters
%     root.peakrate = (rand(1)*14)+6;
%     root.ctrX = (rand(1)*100)+25;
%     root.ctrY = (rand(1)*100)+25;
%     root.sigmaX = 10;
%     root.sigmaY = 10;
%     root.size = 150;
%     root.bins = 50;
%     root.P = P;
%     
%     % simulate ratemap
%     [map] = simulate_ratemap(root);
%     % figure; imagesc(map.z);
% 
%     % generate spiketimes
%     [sim] = simulate_place(map, P(:,1));
%     ST = sim.ST;

    % parameters for egocentric cell
    param.position = P;
    param.ref_point = [75,75];
    param.theta = thetaOptions(i);
    
    [sim] = simulate_ego_cell(param);
    ST = sim.spiketimes;

    % run the model
    [out] = modelMe_V2(P, ST);
    
    % save some stuff for comparison
    g_M(i) = out.model.fitParams.g;
    theta_M(i) = out.model.fitParams.thetaP;
    x_M(i) = out.model.fitParams.xref;
    y_M(i) = out.model.fitParams.yref;
    
    % save some stuff for comparison
    g_CM(i) = out.CM.fitParams.g;
    theta_CM(i) = out.CM.fitParams.thetaP;
    x_CM(i) = out.CM.fitParams.xref;
    y_CM(i) = out.CM.fitParams.yref;
    
    
%     % plot some cosines
%     if mod(i,5) == 0
%         figure(i); hold on;
%         plot(reshape(out.model.Rxyh(5,5,:),10,1));
%         plot(reshape(out.CM.Rxyh(5,5,:),10,1));
%         legend('M', 'CM')
%     end
end

% plot things
figure; hold on;
title('\theta_p: linear v. circular mean');
xlabel('\theta_p (deg)')
histogram(theta_M, linspace(-180,180,10));
histogram(theta_CM, linspace(-180,180,10));
legend('linear', 'circular')

figure; hold on;
title('g: linear v. circular mean');
xlabel('g')
histogram(g_M);
histogram(g_CM);
legend('linear', 'circular')

figure; hold on;
title('y_ref: linear v. circular mean');
xlabel('y_ref')
histogram(y_M);
histogram(y_CM);
legend('linear', 'circular')



