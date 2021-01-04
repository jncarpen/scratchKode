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



