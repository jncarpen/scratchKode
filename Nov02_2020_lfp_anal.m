% what session do you wanna look at?
sessNum = 57;

% grab eeg voltage from a session
eeg_voltage = rawEEG{1,sessNum}(:,1);

% filter the EEG (around the noise frequencies)
nyquistLim = 125; %250/2
inputsig = eeg_voltage;
order = 3;
fcutlow=45; % hz
fcuthigh=55;

[b,a]=butter(order,[fcutlow,fcuthigh]/nyquistLim,'bandpass');
filtsig=filter(b,a,inputsig);  %filtered signal

% find phase angles for each EEG sample
phases = angle(hilbert(filtsig));


for unitNum = 1:length(SpikeTimes{1,sessNum})
    
    % get spiketimes for this unit
    spiketimes = SpikeTimes{1,sessNum}{1,unitNum};

    % grab phase angles for each spikes
    phase_angles = phases(knnsearch(eeg_timestamps, spiketimes));

    % plot a histogram
    n(:, unitNum) = hist(phase_angles, 30);

    % circ_r test (rayleighs test)
    [pval(unitNum), z(unitNum)] = circ_rtest(phase_angles);
end

figure
hist(pval, 30)

figure




% 
imagesc(zscore(n))
% make a histogram of all the p-values (from the circ-r test), and they
% should be uniform (when you pool all of the cells).


% make the periodogram (power spectral density plot)
fig = figure;
set(gcf,'color','w');
periodogram(rawEEG{1,sessNum}{1,1}(:,1),[],5120,250,'power','reassigned')
ax = gca;

%% make reassigned periodogram power spectrum estimate for each session

for sessNum = 61:length(rawEEG)
    % get session type
    ST = sessType{1,sessNum};
    if isempty(ST)
        ST = 'UK';
    end
    
    fig = figure;
    set(gcf,'color','w');
    periodogram(rawEEG{1,sessNum}{1,1}(:,1),[],5120,250,'power','reassigned')
    title(strcat('Session: ', sprintf('%.f', sessNum),', Type: ', ST))
    
    filename = strcat('D:\egoAnalysis\PSD\', sprintf('%.f', sessNum), '_', ST, '.png');
    saveas(fig, filename);
end













