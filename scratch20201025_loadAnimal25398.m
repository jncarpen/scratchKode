%% load animal 25398

recFolderPath = "Z:\jsblacks\Data\25398";
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

% make cell arrays for information you want to pull out
animal_tfilelist = cell(1,length(dirList));
animal_EEGfilelist = cell(1,length(dirList));
animal_posfilelist = cell(1,length(dirList));
animal_trackfilelist = cell(1,length(dirList));
animal_fmEvents = cell(1,length(dirList));

for session= 1:length(dirList)
    folderPath = char(dirList(session));
    
    % grab files with specific extensions
    fileList_t64 = dir(fullfile(folderPath, '*.t64'));
    tfilelist = cell(1,length(fileList_t64));

    fileList_pos = dir(fullfile(folderPath, '*.pos'));
    fileList_tracker = dir(fullfile(folderPath, '*_tracker.mat'));
    % this line was changed from the last animal
    fileList_fmEvents = dir(fullfile(folderPath, 'fmEventsJdp.mat')); 
    
%     % if no file is found...
%     if isempty(fileList_tracker)
%         fileList_tracker = [];
%     end
%     
%     if isempty(fileList_fmEvents)
%         fileList_fmEvents = [];
%     end
    
    % EEG paths
    pathEEG = dir(fullfile(folderPath, '*.eeg'));
    pathEEG2 = dir(fullfile(folderPath, '*.eeg2'));
    pathEEG3 = dir(fullfile(folderPath, '*.eeg3'));
    pathEEG4 = dir(fullfile(folderPath, '*.eeg4'));
    fileList_EEG = [pathEEG, pathEEG2, pathEEG3, pathEEG4];
    EEGfilelist = cell(1,length(fileList_EEG));
        
    for i=1:length(fileList_t64)
        tfilelist{i} = append(fileList_t64(i).folder, '\', fileList_t64(i).name);
    end
    
    for i=1:length(fileList_EEG)
        EEGfilelist{i} = append(fileList_EEG(i).folder, '\', fileList_EEG(i).name); 
    end
    
    for i=1:length(fileList_pos)
        posfilelist{i} = append(fileList_pos(i).folder, '\', fileList_pos(i).name); 
    end
    
     trackerfilelist{1}=[];
    for i=1:length(fileList_tracker)
        if ~isempty(fileList_tracker(i))
            trackerfilelist{i} = append(fileList_tracker(i).folder, '\', fileList_tracker(i).name); 
        else
            trackerfilelist{i} = 'EMPTY';
        end
    end
    
    fmEventsfilelist{1}="OFS";
    for i=1:length(fileList_fmEvents)
        if ~isempty(fileList_fmEvents(i))
            fmEventsfilelist{i} = append(fileList_fmEvents(i).folder, '\', fileList_fmEvents(i).name); 
        else
            fmEventsfilelist{i} = 'OFS';
        end
    end
    
 animal_tfilelist{session} = tfilelist;
 animal_EEGfilelist{session} = EEGfilelist;
 animal_posfilelist{session} = posfilelist;
 animal_trackfilelist{session} = trackerfilelist;
 animal_fmEvents{session} = fmEventsfilelist;
 
%  clear fileList_fmEvents fileList_tracker
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the dataset and save in a matfile

% specify storage location (should be the same for all animals)
storeLoc = "D:\Data\Dataset";
filename = '25398.mat';

% Fill this section out with relevant information for the animal of interest 
recFolderPath = 'Z:\jsblacks\Data\25398'; % change this for each animal (or make it into a function)
LabNotes_flnm = "D:\Data\Labnotes\25398_Oct25.xlsx";

% get length of sessions
totalSessions = length(animal_tfilelist); % grab # of sessions
totalEEGChannels = length(animal_EEGfilelist{1,1}); % grab # EEG channels

% create cell arrays, where each cell is 1 session
spikeTimes = cell(1,totalSessions); % sampling rate of spikeTimes?
unitID = cell(1,totalSessions);
eeg_and_time = cell(1,totalEEGChannels); % [time EEG]
rawEEG = cell(1,totalSessions);
pos = cell(1,totalSessions);
speed = cell(1,totalSessions);
accel = cell(1,totalSessions);
hd = cell(1,totalSessions);
sessInfo = cell(1, totalSessions);
fmEvents = cell(1,totalSessions);


% loop through the filelists for each session
for sessNum = 1:totalSessions
    
    % grab session information (from .pos file) --> getSessInfo.m
     posFilePath = animal_posfilelist{1,sessNum}{1,1};
     [info] = getSessInfo(posFilePath);
     sessInfo{1,sessNum} = info;
    
    % read in files using path locations specified in the animal_xlists
    [spikeTimes{1, sessNum}, unitID{1,sessNum}] = LoadSpikes(animal_tfilelist{1,sessNum});
    eegSession = animal_EEGfilelist{1,sessNum};
    pullPos = io.axona.getPos(animal_posfilelist{1,sessNum}{1,1}); % from BNT* /interpolate NaNs?
    [xInt1, yInt1] = general.interpolatePositions(pullPos(:,1), [pullPos(:,2), pullPos(:,3)]); % interpolate positions for LED1
    [xInt2, yInt2] = general.interpolatePositions(pullPos(:,1), [pullPos(:,4), pullPos(:,5)]); % interp positions for LED2
    intPos = [pullPos(:,1), xInt1, yInt1, xInt2, yInt2]; % squash it back together
    fixedLEDpos = general.fixLedAssignment(intPos); % Fix point assignment to LEDs in tracked positions (computationally heavy)
    
    
    % need to convert to cm (fix this)
    smoothPos(1:length(fixedLEDpos),1) = fixedLEDpos(:,1); % we don't want to smooth time vector
    for col = 2:5
        smoothPos(1:length(fixedLEDpos),col) = general.smoothGauss(fixedLEDpos(:,col), 5); % smooth positions 
        % n.b. sigma is the standard deviation for Gaussian kernel, measured in number of samples (0 = no smoothing).
    end
    pos{1,sessNum} = smoothPos; % save position data for the session
    
    for eegChan = 1:totalEEGChannels
       eeg_and_time{1,eegChan} = read_eeg_file(eegSession{1,eegChan}); % read in all (4) EEG channels
    end
    
   rawEEG{1,sessNum} = eeg_and_time; %% [eeg time]; units? --> match with timestamps *
   
  %% speed + accel + hd
  numLeds = 2; 
  % parse and medfilt pos vector to eliminate major outliers
  t = smoothPos(:,1); % in seconds
  x = medfilt1(smoothPos(:,2)); 
  x2 = medfilt1(smoothPos(:,4));
  y = medfilt1(smoothPos(:,3)); 
  y2 = medfilt1(smoothPos(:,5));
  v = zeros(size(smoothPos,1), numLeds); % velocity
  a = zeros(size(smoothPos,1), numLeds); % acceleration 
 
  for i = 2:size(smoothPos,1)-1
    v(i, 1) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
    v(i, 2) = sqrt((x2(i+1) - x2(i-1))^2 + (y2(i+1) - y2(i-1))^2) / (t(i+1) - t(i-1));
  end
  
  % pad the vector
  v(1,1) = v(2,1);
  v(1,2) = v(2,2);
  v(end, 2) = v(end-1, 2);
  v(end, 1) = v(end-1, 1);

  v1 = v(:,1);
  v2 = v(:,2);
  
  for i = 2:size(v,1)-1
      a(i,1) = (v1(i+1) - v1(i-1))/(t(i+1) - t(i-1));
      a(i,2) = (v2(i+1) - v2(i-1))/(t(i+1) - t(i-1));
  end
         
  speed{1,sessNum} = v; 
  accel{1,sessNum} = a;
  hd{1,sessNum} = rem(atan2d(y2-y, x2-x) + 180, 360); % HD in DEGREES
  
  %% For stuff JS already generated 
  
  % for fmEvents.mat
  if animal_fmEvents{1,sessNum}{1,1} == "OFS"
        fmEvents{1,sessNum} = 'OFS';
  else
      eventStruct = load(animal_fmEvents{1,sessNum}{1,1});
      fmEvents{1,sessNum} = eventStruct;
  end
  
  
end % end of sessions loop





% here's where the code stopped (Oct 25 @ 5:30 PM)
%% pull lab notes 
% modified from the 'pullLabNotes.m' function
% [unitInfo] = pullLabNotes(LabNotes_flnm, 8);

totalTetrodes = 7; % for animal 25398
% there is no tetrode #4 for this animal (not sure why?)
tetNum = ["TT01", "TT02", "TT03", "TT05", "TT06", "TT07", "TT08"]; % add as many as needed...

if totalTetrodes > length(tetNum)
    disp("Error: Add options for additional tetrodes in getUniqueID.m")
end

unitInfo = []; % initialize empty 'table'

for tetSheet = 1:totalTetrodes
tetInfo = readtable(LabNotes_flnm, 'Sheet', tetNum(tetSheet));

    cRows = [];
    for row = 1:height(tetInfo)
        if isnan(tetInfo.Cluster(row))
            cRows(row) = row;
        end
    end
    cRows = nonzeros(cRows); 
    tetInfo(cRows,:) = []; % delete 'c' rows
    
    % Columns to keep 
    tetInfoTite = table(tetInfo.CutFile, tetInfo.TypeOfSession, tetInfo.Tetrode, tetInfo.Cluster, tetInfo.UniqueID, tetInfo.ID);
    tetInfoTite.Properties.VariableNames = {'CutFile', 'TypeOfSession', 'Tetrode', 'Cluster', 'UniqueID', 'ID'};

    for row = 1:height(tetInfoTite)
        cutFl = tetInfoTite.CutFile{row,1}; % grab cutfile for current row
        cutFl = extractAfter(cutFl,"25398\"); % I want to remove the .clusters at the end but not every string has it
        if strfind(cutFl,".clusters") == 13
            cutFl = extractBefore(cutFl,".clusters");
        end
        
        if strfind(cutFl, "_") == 11 % in position 11?
            cutFl = extractBefore(cutFl, "_");
        end
        tetInfoTite.CutFile{row,1} = cutFl;
    end
    
    unitInfo = [unitInfo; tetInfoTite]; % concat tables from multiple tetrodes together
end

%% get unique id #s for this animal
  [uniqueID, neuronID, sessType] = getUniqueID(unitInfo, animal_tfilelist);

%% Bin stuff
spikeTrain = cell(1, totalSessions);
for sessNum = 1:totalSessions
    if ~isempty(pos{1,sessNum}) && ~isempty(spikeTimes{1,sessNum}) % do we have everything we need?
        spikeTrain{1,sessNum} = binSpikes(pos{1,sessNum}(:,1), spikeTimes{1,sessNum}); 
    else
        spikeTrain{1,sessNum} = []; % no tracker file/spikeTimes for that session...
    end
end
    
    %% Remove cells that don't have a unique cell ID 
    SpikeTrain = cell(1,totalSessions); % capital S
    SpikeTimes = cell(1,totalSessions);
    UniqueID = cell(1,totalSessions);
    UnitID = cell(1,totalSessions);
    
    for sessNum = 1:totalSessions
        STRN = spikeTrain{1,sessNum};
        STM = spikeTimes{1,sessNum};
        U = uniqueID{1,sessNum};
        ID = unitID{1,sessNum};
        keep = find(~cellfun(@isempty,U)); % find cells w/ unique IDs
        SpikeTrain{1,sessNum}  = STRN(keep);
        SpikeTimes{1,sessNum} = STM(keep);
        UniqueID{1,sessNum} = U(keep);
        UnitID{1,sessNum} = ID(keep);
    end
    

%% SAVE DATASET
% save D:\Data\Dataset\25398_v1.mat SpikeTimes SpikeTrain UnitID rawEEG pos speed accel hd sessInfo fmEvents unitInfo UniqueID neuronID sessType
save('D:\Data\Dataset\25398_eeg.mat', 'rawEEG', '-v7.3')


%% CORRECT POSITION
% now im going to correct the position file

pos_cm = cell(1,length(pos));
for sessNum = 1:length(pos)
    pos_cm{1,sessNum} = Position_cm(pos, sessInfo, sessNum);
end



% save the dataset again
save D:\Data\Dataset\25398\25398_v2.mat SpikeTimes SpikeTimes_thresh SpikeTrain UnitID pos pos_cm speed speed_cm accel hd sessInfo fmEvents unitInfo UniqueID neuronID sessType hwCoord boxCtr
