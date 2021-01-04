% looking at session # 57 (fm(13:16,:)

% what session do you want to look at?
sessNum = 57;

% grab timestamps for session
t = pos_cm{1,sessNum}(:,1);
fs = mode(diff(t));

% grab foster maze events for this session
events = fmEvents{1,sessNum}.fmEvents;
times = events.times;
type = events.type;

% grab EEG for this session
chanNum = 3;
eeg_val = rawEEG{1,sessNum}{1,chanNum}(:,1);
eeg_time = rawEEG{1,sessNum}{1,chanNum}(:,2);
eeg_fs = mode(diff(eeg_time));

% find timestamps for each home drinking event
type_count = count(type, "DRINKING_HOME");
home = times.*type_count; % pull out drinking_random events
home = home(any(home,2),:); % pulls out 0 values
home_ts = knnsearch(t, home(:,1));
home = t(knnsearch(t, home(:,1))); % match with closest timestamp

% find decision point
type_count = count(type, "DECISION_POINT");
decision = times.*type_count; % pull out drinking_random events
decision = decision(any(decision,2),:); % pulls out 0 values
decision_ts = knnsearch(t, decision(:,1)); 
decision = t(knnsearch(t, decision(:,1))); % match with closest timestamp

% grab decision point (only first column)
type_count = count(type, "DRINKING_RANDOM");
random = times.*type_count; % pull out drinking_random events
random = random(any(random,2),:); % pulls out 0 values
random = knnsearch(t, random(:,1)); 
random = t(knnsearch(t, random(:,1))); % match with closest timestamp


% match behavioral timestamps with eeg timestamps
match_eeg_timestamp = knnsearch(eeg_time,  decision(:,1));

numSamples = 5000; %(+/- 2 seconds)
numSeconds = eeg_fs * (numSamples/2);
figure
set(gcf,'color','w');
for eventNum = 1:length(match_eeg_timestamp)
    ts = match_eeg_timestamp(eventNum);
    plot((-numSeconds:eeg_fs:numSeconds).*1000, 100*eventNum + eeg_val(ts:ts+numSamples), 'LineWidth', .25);
    hold on
end

xline(0, 'LineWidth', 1.5)
xlim([-numSeconds*1000 numSeconds*1000])
xlabel('time (ms)')
box off
yticks([])

for i =1:length(random)
xline(random(i), 'color', 'green');
hold on
end

for i =1:length(decision)
xline(random(i), 'color', 'red');
hold on
end

for i =1:length(home)
xline(random(i), 'color', 'cyan');
hold on
end

%% visualize EEG (with events)


numChan = length(rawEEG{1,sessNum}); % grab # of channels

    for ts = 1:1000:length(rawEEG{1,sessNum}{1,1})
        for i=1:numChan

            plot(rawEEG{1,sessNum}{1,1}(ts:ts+1000,2),100*i+rawEEG{1,sessNum}{1,i}(ts:ts+1000,1),'k')
            hold on
            title("vizEEG")
            xlabel("time (s)")
            ylabel("amplitude")
            
        end
        
        events_now = home(find(home>rawEEG{1,sessNum}{1,1}(ts,2) & home>rawEEG{1,sessNum}{1,1}(ts+1000,2)));
        if ~isempty(events_now)
            for i =1:length(events_now)
            xline(events_now(i), 'color', 'red');
            hold on
            end
        end
            
        pause
        clf
        
    end
