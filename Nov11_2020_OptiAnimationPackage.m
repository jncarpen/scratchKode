
% load data
load data.mat

% choose a session
sessNum = 1;

% display session type
disp(strcat("session type: ", convertCharsToStrings(sessType{1,sessNum})))

% get csv filepath
csv_filepath = optitrack_directory{1,sessNum};


for sessNum = 1:length(optitrack_directory)
    csv_filepath = optitrack_directory{1,sessNum};
    if ~isempty(csv_filepath)
        csv_filepath = convertCharsToStrings(csv_filepath);
        copyfile (csv_filepath, 'D:\Data\Dataset\25398\OptiAnimation_Package\CSV')
    end
end