function spike_log_out = clean_spike_map(spike_log)

% Remove neurons that did not meet criteria post-examination
spike_log_out = spike_log(spike_log.useNeuron == 1,:);

end