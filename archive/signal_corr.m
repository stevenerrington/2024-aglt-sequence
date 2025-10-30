session_list = unique(spike_log.session);
count = 0;
rng(20)

areal_corr = [];
for session_i = 1:length(session_list)

    session_entries = []; session_entries = find(strcmp(session_list{1},spike_log.session));
    session_aud_neurons = []; session_aud_neurons = find(strcmp(spike_log.area(session_entries,:),'auditory'));
    session_frontal_neurons = []; session_frontal_neurons = find(strcmp(spike_log.area(session_entries,:),'frontal'));

    [dual_neuron_n, min_idx] = min([length(session_aud_neurons), length(session_frontal_neurons)]);

    switch min_idx
        case 1
            session_frontal_neurons = randsample(session_frontal_neurons, dual_neuron_n);
        case 2
            session_aud_neurons = randsample(session_aud_neurons, dual_neuron_n);
    end

    % correlation
    for pair_i = 1:dual_neuron_n
        count = count + 1;
        aud_sdf = []; aud_sdf = pca_sdf_out(session_aud_neurons(pair_i),1000+[0:2650]);
        frontal_sdf = []; frontal_sdf = pca_sdf_out(session_frontal_neurons(pair_i),1000+[0:2650]);

        [xcorr_out(count,:), lag] = xcorr(aud_sdf, frontal_sdf,'normalized');

    end

end

figuren;
plot(lag,mean(xcorr_out))
xlim([-500 500])