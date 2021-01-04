% remove duplicate fmEvents    

% grab a random session
kk = 3;
type_raw = fmEvents{1,kk}.fmEvents.type; % example cell array

type = type_raw;
for i = 1:length(type_raw)-1
    if isequal(type_raw{i,1}, type_raw{i+1,1})
        type{i+1,1} = [];
    end  
end

% get rid of empty cells
index = cellfun(@isempty, type) == 0;
type = type(index);

% now this can be incorporated back into the getHomeRuns.m script