%% RUN MODEL ON OPTITRACK 2D DATA
% Jan 24, 2020 (Sunday)
% Jordan Carpenter.
% LEGEND -
% COL 1:        frame # (integer)- starts at 0
% COL 2:        time (seconds)
% COL 3:        X-rotation (quaternions)
% COL 4:        Y-rotation (quaternions)
% COL 5:        Z-rotation (quaternions)
% COL 6:        W-rotation (quaternions)
% COL 7:        X-position (meters)
% COL 8:        Y-position (meters)
% COL 9:        Z-rotation (meters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in the dataset
load('D:\Data\Dataset\opti4JS\optiPackJan182021.mat')
totalSessions = length(optiTrack);



for sess = 1:totalSessions
    t = optiTrack{1,sess};
    if ~isempty(table)
        clear P ST HD
        t = optiTrack{1,sess}; isempty(t)
        % get variables from csv & convert from Y-Up RH system
        self.timestamp = table2array(t(:,2));
        self.q.x = table2array(t(:,3));
        self.q.y = table2array(t(:,5));
        self.q.z = table2array(t(:,4));
        self.q.w = table2array(t(:,6));
        self.position.x = table2array(t(:,7));
        self.position.y = table2array(t(:,9));
        self.position.z = table2array(t(:,8));
        self.fs = round(1/mode(diff(self.timestamp))); % in hz
        self.n_frame = numel(self.timestamp);
        
        % make x-values positive
        self.q.x = -self.q.x;
        self.position.x = -self.position.x;
        
        % calculate movement speed, direction, and angle
        self.movement = calculate_movement(self);
        
        % calculate rotation matrix
        self.rotation_matrix = calculate_rotation_matrix(self.q);
        
        % calculate yaw (azimuth), pitch and roll allocentric head direcitons:
        self.direction = calculate_direction(self);
        
        % head direction (yaw in radians, -pi:pi)
        HD = rad2deg(self.direction.azimuth);
        
        % position vector (meters)
        P = [self.timestamp, self.position.x, self.position.y];
        
        % simulate a cell
        param.position = P;
        param.ref_point = [0,0];
        param.theta = 270;
        param.hd = HD;
        param.ctr_mass = [0,0];
        param.noise = 0;
        param.width = 10;
        
%         [sim] = simulate_place_egoMod(param);

        [sim] = simulate_ego_cell(param);
        ST = sim.spiketimes;
        pathPlot_hd(P, ST, HD)
%         [sim] = simulate_place_egoMod(param);
        
        [out] = modelMe_V2(P, ST, HD);
        plot_vectorMod(out)
        figure; imagesc(flipud(out.data.rxyS));
        
        
        
        
        
        
        
        % find number of units in the session
        howManyUnits = length(STV{1,sess});
        for unit = 1:howManyUnits
            % grab spiketimes 
            ST = STV{1,sess}{1,unit};
            
            % @criteria: unit must spike more than 150x
            if length(ST) > 150
                
                % run the model
                [out] = modelMe_V2(P, ST, HD);

            end
        end
    end
end


% look at stuff
g = out.model.fitParams.g;
thetaP = deg2rad(out.model.fitParams.thetaP);
xref = out.model.fitParams.xref;
yref = out.model.fitParams.yref;
bins = linspace(-180,180,10);
avgCos = 1 + g*cos(bins-thetaP);

figure; hold on;
plot(bins, avgCos);

clear M D
for i = 1:10
    M(i) = nanmean(out.model.Rxyh(:,:,i), 'all');
    D(i) = nanmean(out.data.Rxyh(:,:,i), 'all');
end

figure; hold on;
plot(bins, M);
plot(bins,D, '.')















