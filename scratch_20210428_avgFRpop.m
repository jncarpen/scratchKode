pop = place1000;

for i = 1:length(pop)
    t = pop(i).param.P(:,1); % timestamps
    tpf = mode(diff(t)); % time per frame
    ST = pop(i).ST; % spiketimes
    [spikeTrain, ~] = binSpikes(t, ST); % spiketrain (counts)
    inst_fr = imgaussfilt(spikeTrain./tpf,50); 
    avg_fr(i) = mean(inst_fr);
    max_fr(i) = max(inst_fr);
end

%% detect place fields

ctrMassX = zeros(length(pop),1)*nan;
ctrMassY = zeros(length(pop),1)*nan;
mean_infield_rate = zeros(length(pop),1)*nan;
totalArea = zeros(length(pop),1)*nan;
for i = 1:length(pop)
    P = pop(i).param.P;
    ST = pop(i).ST;
    map = analyses.map(P,ST);
    [fieldMap, fields] = analyses.placefield(map, 'threshold', 0.2, 'minBins', 20, 'minPeak', 0.3);
    borderCov(i) = analyses.borderCoverage(fields);
    % find the largest field
    if ~isempty(fields)
        for fieldNow = 1:length(fields)
            fieldArea(fieldNow) = fields(fieldNow).area;
        end
        totalArea(i) = sum(fieldArea);
        [~, bigField] = nanmax(fieldArea);
        ctrMassX(i) = fields(bigField).x.*2.5;
        ctrMassY(i) = fields(bigField).y.*2.5;
        mean_infield_rate(i) = fields(bigField).meanRate;
        clear fieldArea
    end
end
% 
% 
% 
% 
% %% distance from nearest boundary
% % x = ctrMassX; y = ctrMassY;
% % Q1 = find(x>=75 & x<150 & y>=75 & y<150);
% % Q2 = find(x>=75 & x<150 & y>=0 & y<75);
% % Q3 = find(x>=0 & x<75 & y>=0 & y<75);
% % Q4 = find(x>=0 & x<75 & y>=75 & y<150);
% 
% 
% 















