% november 5, 2020
% description: plot the trajectory for one trial, color-coded by task-phase.

% choose a session # 
sessNum = 57; % fm14

% get position information
P = pos_cm{1, sessNum};
x = P(:,2); y = P(:,3);
t = P(:,1); fs = mode(diff(t));

% get hw location
x_ref = hwCoord{1,sessNum}(1,1);
y_ref = hwCoord{1,sessNum}(1,2);

% get task-phase information
events = fmEvents{1,sessNum}.fmEvents; % struct
type = events.type;
time_start = events.times(:,1);
time_stop = events.times(:,2);

% find event
phase1_idx = find(count(type, "DRINKING_HOME") == 1);
% n = floor(length(phase1_idx)/2);
n = 6;
phase1_idx = phase1_idx(n);

if type(phase1_idx) == "DRINKING_HOME" & type(phase1_idx+1) == "DRINKING_RANDOM" & type(phase1_idx+2) == "DECISION_POINT" & type(phase1_idx+3) == "DRINKING_HOME"
    
    % get info for trial of interest
    home_times = t(knnsearch(t, time_start(phase1_idx))):fs:t(knnsearch(t, time_stop(phase1_idx)));
    idx_home = knnsearch(t, home_times');
    
    rand_times = t(knnsearch(t, time_stop(phase1_idx))):fs:t(knnsearch(t, time_start(phase1_idx+2)));
    idx_rand = knnsearch(t, rand_times');
    
    dec_times = t(knnsearch(t, time_start(phase1_idx+2)));
    idx_dec = knnsearch(t, dec_times');
    
    home2_times = t(knnsearch(t, time_start(phase1_idx+2)))+fs:fs:t(knnsearch(t, time_stop(phase1_idx+3)));
    idx_home2 = knnsearch(t, home2_times');
    
    % plot
    figure(1)
    set(gcf,'color','w');
    hold on;
    plot(x(idx_home), y(idx_home), 'LineWidth', 1.25, 'color', [1 .75 .5])
    plot(x(idx_rand), y(idx_rand), 'LineWidth', 1.25, 'color', [0 .5 1])
    plot(x(idx_home2), y(idx_home2), 'LineWidth', 1.25, 'color', [.8 0 1])
    
    % add decision point
    d_pnt = plot(x(idx_dec), y(idx_dec),  '*', 'MarkerSize', 25);
    set(d_pnt, 'markerfacecolor', 'blue');
    
    % add home well location
    hw_loc = plot(x_ref, y_ref, 'o', 'MarkerSize', 8);
    set(hw_loc, 'markerfacecolor', 'red');
    
    % random well location
    rand_loc = plot(x(idx_rand(end)), y(idx_rand(end)), 'o', 'MarkerSize', 8);
    set(rand_loc, 'markerfacecolor', 'blue');
    
    % add legend
    l = legend('home drinking', 'random foraging', 'return home 2', 'decision pt', 'home well', 'random well')
    l.FontName = 'Calibri Light'; l.FontSize = 14; l.Location = 'northeastoutside';
    
    xlim([0 150]); ylim([0 150])
    xlabel('x (cm)', 'FontName', 'Calibri Light', 'FontSize', 20);
    ylabel('y (cm)', 'FontName', 'Calibri Light', 'FontSize', 20);
    
    title('A25398, Foster Maze + Object (S57)', 'FontName', 'Calibri Light', 'FontSize', 40, 'FontWeight', 'Normal')
    
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Calibri Light','fontsize',15)
    b = get(gca,'YTickLabel');
    set(gca,'YTickLabel',b,'FontName','Calibri Light','fontsize',15)

else
    disp("Error: change value of 'n'")
end











%% SCRATCH

% home_start = t(knnsearch(t, time_start(find(count(type, "DRINKING_HOME") == 1))));
% home_stop = t(knnsearch(t, time_stop(find(count(type, "DRINKING_HOME") == 1))));
% 
% random_start = t(knnsearch(t, time_start(find(count(type, "DRINKING_RANDOM") == 1))));
% random_stop = t(knnsearch(t, time_stop(find(count(type, "DRINKING_RANDOM") == 1))));
% 
% decision = t(knnsearch(t, time_start(find(count(type, "DECISION_POINT") == 1))));
% 
% % plot the nth trial
% n = floor(length(decision)/2); % find somewhere in the middle
% 
% figure(1)
% hold on;
% plot(x(find(t == home_start(n)):find(t == home_stop(n))), y(find(t == home_start(n)):find(t == home_stop(n))), 'LineWidth', 1.5, 'color', 'red')
% plot(x(find(t == random_start(n)):find(t == random_stop(n))), y(find(t == random_start(n)):find(t == random_stop(n))), 'LineWidth', 1.5, 'color', 'blue')
