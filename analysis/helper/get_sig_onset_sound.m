function sig_onset_times = get_sig_onset_sound(ops, sound_sdf, neuron_baseline)

sound_sdf_window = ops.sound_sdf_window;
sound_sdf_average = nanmean(cell2mat(sound_sdf(:,1)));
baseline_mean = nanmean(neuron_baseline);
baseline_std = nanstd(neuron_baseline);

sig_idx = sound_sdf_window(sound_sdf_average > baseline_mean+(baseline_std*ops.sig_threshold) |...
    sound_sdf_average < baseline_mean-(baseline_std*ops.sig_threshold));


[start, len] = ZeroOnesCount(ismember(ops.sound_sdf_window,sig_idx));
idx_w_sigTimes = find(len > ops.min_sig_time);

sig_onset_times = sound_sdf_window(start(idx_w_sigTimes));


end