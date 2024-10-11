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
onset_times = [0	563	1126	1689	2252];
offset_times = [413	976	1539	2102	2665];
sound_label = {'YAG'	'LEK'	'PAV'	'RAP'	'KEM'};

sound_i = 14;

audio_file = stimuli_table.filename{file_i};

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
    box off
    set(gca,'XColor',[0 0 0], 'YColor', [1 1 1])
    title(sound_label{search_win_idx})
    xlim([0 400])
end



