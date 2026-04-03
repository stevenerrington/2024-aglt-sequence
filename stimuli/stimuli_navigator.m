stimuli_table = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t_stimuli'));

n_stimuli = size(stimuli_table,1);
figuren('Renderer', 'painters', 'Position', [100 100 400 1200]); hold on;;

stimuli_onset_search_start = [0.0001 0.55 1.1 1.65 2.2];

threshold = 0.05;

for file_i = 1:n_stimuli

    audio_file = stimuli_table.filename{file_i};

    [sounddata,fs_check] = audioread...
        (['C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\stimuli\' audio_file]);

    sounddata = sounddata(:,1);


    for search_win_idx = 1:5
        search_win = round([stimuli_onset_search_start(search_win_idx),...
            stimuli_onset_search_start(search_win_idx)+0.25]*fs_check);
        sounddata_snippet = sounddata(search_win(1):search_win(2));

        startingIndex(search_win_idx) = find(abs(sounddata_snippet > threshold),1,'first');

        sound_time(search_win_idx,1) = (startingIndex(search_win_idx)/...
            fs_check)+stimuli_onset_search_start(search_win_idx);

    end


    time = (1:length(sounddata))/fs_check;

    clip_duration(file_i,1) = max(time);

    subplot(n_stimuli,1,file_i)
    plot(time, sounddata(:,1))
    vline(sound_time, 'r-')


    sound_time_out(file_i,:) = sound_time*1000';
end


%%
onset_times = [0	563	1126	1688	2252];
offset_times = [413	976	1539	2102	2665];
sound_label = {'YAG'	'LEK'	'PAV'	'RAP'	'KEM'};

sound_i = 14;

audio_file = stimuli_table.filename{file_i};

[sounddata,fs_check] = audioread...
    (['/Users/stevenerrington/Desktop/Projects/2024-aglt-sequence/stimuli/' audio_file]);

sounddata = sounddata(:,1);

time =((1:length(sounddata))/fs_check)*1000;

colorscale = cbrewer('qual','Set2',5);

figuren('Renderer', 'painters', 'Position', [100 100 1200 200]); hold on;

count = 0;
sound_out = [];

for search_win_idx = [1 2 5 3 4]
    count = count + 1;
    subplot(1,5,count)
    idx = []; idx = find(time > onset_times(search_win_idx) & time < offset_times(search_win_idx));

    idx = idx(1:9106);
    plot(([1:length(idx)]/fs_check)*1000, sounddata(idx),'color',colorscale(search_win_idx,:));

    sound_out(search_win_idx,:) = sounddata(idx);

    box off
    set(gca,'XColor',[0 0 0], 'YColor', [1 1 1])
    title(sound_label{search_win_idx})
    xlim([0 400])
end


%%
x = sound_out(search_win_idx,:);
fs = fs_check;

% Parameters (adjust if desired)
winLength = 1024;        % Window size (samples)
hopLength = winLength/2; % Hop size (overlap = 50%)
nfft = 1024;

window = hamming(winLength, 'periodic');

numFrames = floor((length(x) - winLength) / hopLength) + 1;
S = zeros(nfft/2 + 1, numFrames);

% Short-Time Fourier Transform
for k = 1:numFrames
    idx = (k-1)*hopLength + (1:winLength);
    frame = x(idx) .* window;

    X = fft(frame, nfft);
    S(:, k) = abs(X(1:nfft/2 + 1));
end

% Axes
timeAxis = (0:numFrames-1) * hopLength / fs;
freqAxis = (0:nfft/2) * fs / nfft;

% Plot
figure;
imagesc(timeAxis, freqAxis, 20*log10(S + eps));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram');
colorbar;