%% February 15, 2021
% Description: Load Adrien Peyrache's TH-1 dataset (spikes + position)

%%
% foldername for a single session
folderPath = 'D:\Data\Dataset\sample data\Peyrache\Mouse12-120806.tar\Mouse12-120806\Mouse12-120806';
TH1data.animal = 'Mouse12';
TH1data.session = '120806';

clear spk_directory posfile_path
session_directory = dir(fullfile(folderPath));
for fl = 1:length(session_directory)
    if contains(session_directory(fl).name, '.res')
        cellidx = str2num(extractAfter(session_directory(fl).name,19)); % extract unit number
        spk_directory(cellidx).path = [session_directory(fl).folder, '\', session_directory(fl).name];
    elseif contains(session_directory(fl).name, '.whl')
        % there should only be one .whl file
        posfile_path = [session_directory(fl).folder, '\', session_directory(fl).name];
    end
end

% save the position file
position = load(posfile_path);
position(position == -1) = nan;

% load & save unit spike times (seconds)
for unit_num = 1:length(spk_directory)
    TH1data.unit(unit_num) = 
end


