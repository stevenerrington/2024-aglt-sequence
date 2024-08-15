% Now reference Google sheet
spike_log = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t_spikes'));

count = 0 ;
for session_i = 1:size(ephysLog,1)
    outfile_name = ephysLog.session{session_i}; % Processed file name
    load(fullfile(dirs.mat_data,[outfile_name '.mat']),'spikes');

    sessionidx = []; sessionidx = find(strcmp(spike_log.session,outfile_name));
    session_wav_labels = []; session_wav_labels = spike_log.unitWAV(sessionidx);

    clear waves mean_wave std_wave

    for neuron_i = 1:length(sessionidx)
        count = count + 1;
        waves = spikes.waveform.(session_wav_labels{neuron_i});
        mean_wave = nanmean(waves);
        std_wave = nanstd(nanstd(waves-mean_wave));

        mean_wave_out{count} = waves;
        SNR(count,1) = (max(mean_wave) - min(mean_wave))./(2 * std_wave);
    end

end

figuren;
histogram(SNR,100)

sum(SNR < 10)

SNR < 10