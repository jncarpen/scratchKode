%% February 15, 2021
% Description: Load Adrien Peyrache's TH-1 dataset (position files only)

% make a directory of all position files
folderpath = 'D:\Data\External Data\Peyrache-TH1\PositionFiles';
posfile_dir = dir(fullfile(folderpath));

% load and save each position file
count = 1;
% clear thdata
for fl = 32:length(posfile_dir)
    if contains(posfile_dir(fl).name, 'Mouse')
        pathnow = [posfile_dir(fl).folder, '\', posfile_dir(fl).name];
        position_now = load(pathnow);
        position_now(position_now == -1) = nan; % set -1 to nan
        thdata(count).P = position_now;
        thdata(count).path = pathnow;
        count = count + 1;
    end
end

save('D:\Data\External Data\Peyrache-TH1\PositionFiles\thdata.mat', 'thdata', ...
    '-v7.3');

