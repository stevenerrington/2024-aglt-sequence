
%%
onset_times = [0	563	1126	1689	2252];
offset_times = [413	976	1539	2102	2665];
sound_label = {'YAG'	'LEK'	'PAV'	'RAP'	'KEM'};

sound_i = 14;

audio_file = stimulusLog.filename{sound_i};

[sounddata,fs_check] = audioread...
    (['C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\stimuli\' audio_file]);

sounddata = sounddata(:,1);

time =((1:length(sounddata))/fs_check)*1000;

colorscale = cbrewer('qual','Set2',5);

figuren('Renderer', 'painters', 'Position', [100 100 1200 200]); hold on;

count = 0;

for search_win_idx = [1 2 5 3 4]
    count = count + 1;
    subplot(1,5,count)
    idx = []; idx = find(time > onset_times(search_win_idx) & time < offset_times(search_win_idx));
    plot(([1:length(idx)]/fs_check)*1000, sounddata(idx),'color',colorscale(search_win_idx,:))

    sound_data(search_win_idx,:) = sounddata(idx(1:9106));
    box off
    set(gca,'XColor',[0 0 0], 'YColor', [1 1 1])
    title(sound_label{search_win_idx})
    xlim([0 400])
end

fs = fs_check;


% Spectrogram parameters
windowLength = 128;     % samples
overlap = 8;           % samples
nfft = 256;

% Plot spectrogram
figure('Renderer','painters','Position',[200 200 1500 200]); hold on;;

for search_win_idx = [1 2 5 3 4]
    subplot(1,5,search_win_idx)
    spectrogram(sound_data(search_win_idx,:), windowLength, overlap, nfft, fs, 'yaxis');
    colorbar;
    clim([-100 -30])
    ylim([0 10])
    colormap(slanCM('thermal',200))
    title(sound_label{search_win_idx})
    xlim([0 400])


end