% add path where log files are being stored
addpath(genpath('D:\Idun\log'))

% folder path
filepath = 'D:\Idun\log\20210131';

% directory
log_dir = dir(filepath);

clear dirList
count = 1;
for dirIndex = 1:length(log_dir)
    L = strlength(log_dir(dirIndex).name);
    if L > 2
        dirList{count} = append(log_dir(dirIndex).folder, '\', log_dir(dirIndex).name);
        count = count + 1;
    end
end

clear catFP
for i = 1:length(dirList)
    now = load(dirList{1,i});
    catFP(i).ST = now.FP.ST;
    catFP(i).modelData = now.FP.modelData;
    catFP(i).ve_place = now.FP.ve_place;
    catFP(i).ve_rh = now.FP.ve_rh;
    catFP(i).ts_hd = now.FP.ts_hd;
    catFP(i).ts_rh = now.FP.ts_rh;
    catFP(i).shuffErr = now.FP.shuffErr;
    catFP(i).g = now.FP.g;
    catFP(i).thetaP = now.FP.thetap;
    catFP(i).xref = now.FP.xref;
    catFP(i).yref = now.FP.yref;
    
    dist = catFP(i).ve_rh;
end
