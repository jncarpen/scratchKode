% grab optitrack files

animal = 25398;
recFolderPath = strcat('D:\Data\', sprintf('%.f', animal));
recDir = dir(recFolderPath);
dirList = cell(1, length(recDir));
    
count = 1;
for dirIndex = 1:length(recDir)
    L = strlength(recDir(dirIndex).name);
    % check that folder is named as a session
    if L == 10 
        % squash strings together to get correct filepath
        dirList{count} = append(recDir(dirIndex).folder, '\', recDir(dirIndex).name);
        count = count + 1;
    end
end

dirList = dirList(~cellfun(@isempty, dirList)); % clear empty cells

%% get correct order for dirList (sorted by date)

MDYS = cell(length(dirList), 1);

for fileNum = 1:length(dirList)
    splitName = split(dirList{1,fileNum},'\');
    serial = splitName{end,1}; % grab serial #
    MDYS{fileNum, 1} = strcat(serial(3:4), serial(1:2), serial(5:8), serial(9:10));
end

% sort by mo-day-yr-session
[~, I] = sort(MDYS); % grab order indices
dirList = dirList(I); 

%% create an empty array to store the paths to the optitrack (.csv) files

animal_optilist = cell(1,length(dirList));

for sessNum= 1:length(dirList)
    folderPath = char(dirList{1,sessNum});
    
    % grab files with specific extensions
    fileName_csv = dir(fullfile(folderPath, '*.csv'));
    
    if length(fileName_csv) > 1
        fileName_csv = fileName_csv(1); % take the first (for now)
    end
    
    
    if ~isempty(fileName_csv)
        % save information for this session
        animal_optilist{1, sessNum} = strcat(fileName_csv.folder, '\', fileName_csv.name);
    else
        animal_optilist{1, sessNum} = [];
    end
    
    clear fileList_csv
end


%% now make it into a cell array
totalSessions = length(animal_optilist);
optiTrack = cell(1,totalSessions);

for sessNum = 1:totalSessions
    csv_path = animal_optilist{1,sessNum};
    
    if ~isempty(csv_path)
        optiTrack{1,sessNum} = readtable(csv_path, 'HeaderLines',7); 
    else
        optiTrack{1,sessNum} = [];
    end
end

%% make a struct for each session
clear optiData
for i = 1:length(optiTrack)
    if ~isempty(optiTrack{1,i})
        % grab the table for this session
        clear self
        t = optiTrack{1,i};
        
        % Optitrack Motive use a Y-up right handed coordinate system, which
        % means Y is up and positive X axis is to the left in the original data.
        % I prefer to use a Z-up left handed system, which means Z is up and
        % positive X axis is to the right. So, while we assign the variables we also
        % translate between these reference frames; subtracting the x axis to get a
        % positive x to the right, and switching out the Y and Z axis:
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
        
        self.q.x = -self.q.x;
        self.position.x = -self.position.x;
        
        % calculate movement speed, direction, and angle
        self.movement = calculate_movement(self);
        
        % calculate rotation matrix
        self.rotation_matrix = calculate_rotation_matrix(self.q);
        
        % calculate yaw (azimuth), pitch and roll allocentric head direcitons:
        self.direction = calculate_direction(self);
        
        % simplify calculations by referencing (put this into the self.movement struct?)
        dv_x = squeeze(self.rotation_matrix(:, 1, :))';
        dv_y = squeeze(self.rotation_matrix(:, 2, :))';
        dv_z = squeeze(self.rotation_matrix(:, 3, :))';
        
        x = self.position.x;
        y = self.position.y;
        z = self.position.z;

        % Add position values to unit vectors for animating the rigid body axis at
        % the location of the animal in 3D.
        unit_vec_length = 6;  % length of the unit vector in animation.

        self.x_vector.x = [x, x + dv_x(:, 1) * unit_vec_length];
        self.x_vector.y = [y, y + dv_x(:, 2) * unit_vec_length];
        self.x_vector.z = [z, z + dv_x(:, 3) * unit_vec_length];

        self.y_vector.x = [x, x + dv_y(:, 1) * unit_vec_length];
        self.y_vector.y = [y, y + dv_y(:, 2) * unit_vec_length];
        self.y_vector.z = [z, z + dv_y(:, 3) * unit_vec_length];

        self.z_vector.x = [x, x + dv_z(:, 1) * unit_vec_length];
        self.z_vector.y = [y, y + dv_z(:, 2) * unit_vec_length];
        self.z_vector.z = [z, z + dv_z(:, 3) * unit_vec_length];

        movement_vector = ...
            [cos(self.movement.direction), sin(self.movement.direction)] * 10;

        self.m_vector.x = [x, x + movement_vector(:, 1)];
        self.m_vector.y = [y, y + movement_vector(:, 2)];
        self.m_vector.z = [z, z];

        % save the struct
        optiData{1,i} = self;
        
    else
        optiData{1,i} = [];
    end
end



%% match up the spiketimes
clear ST
for i = 1:length(SpikeTimes)
    if ~isempty(optiTrack{1,i})
        tNow = optiData{1,i}.timestamp;
    else
        tNow = pos{1,i}(:,1);
    end
    
        for j = 1:length(SpikeTimes{1,i})
            stNow = SpikeTimes{1,i}{1,j};
            stAdj = tNow(knnsearch(tNow, stNow));
            ST{1,i}{1,j} = stAdj;
        end
end



%% save stuff
save('D:\Data\Dataset\25398\25398OT.mat', 'optiData', 'ST', '-v7.3')



















