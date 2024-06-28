function signal_out = adjustBaseline(signal_in,ops)

baseline_idx = find(ismember(ops.timewin,ops.baseline));

for ch_i = 1:size(signal_in,1)
    for trial_i = 1:size(signal_in,3)
        signal_out(ch_i,:,trial_i) = signal_in(ch_i,:,trial_i) - mean(signal_in(ch_i,baseline_idx,trial_i));
    end
end
