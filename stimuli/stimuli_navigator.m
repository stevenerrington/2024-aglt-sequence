stimuli_table = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t_stimuli'));

n_stimuli = size(stimuli_table,1);
figuren('Renderer', 'painters', 'Position', [100 100 400 1200]); hold on;;

stimuli_onset_search_start = [0.0001 0.55 1.1 1.65 2.2];

threshold = 0.05;

for file_i = 1:n_stimuli

    audio_file = stimuli_table.filename{file_i};

    [sounddata,fs_check] = audioread...
        (['C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\stimuli\' audio_file]);

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


