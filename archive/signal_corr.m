session_list = unique(spike_log.session);
count = 0;
rng(20)

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

    timewin = [0:413];

    % correlation
    for pair_i = 1:dual_neuron_n
        count = count + 1;
        areal_corr_baseline(count,:) = corr(pca_sdf_out(session_frontal_neurons(pair_i),1000+[-413:0])', pca_sdf_out(session_aud_neurons(pair_i),1000+[-413:0])');
        areal_corr_signal_1(count,:) = corr(pca_sdf_out(session_frontal_neurons(pair_i),1000+[sound_onset_ms(1):sound_onset_ms(1)+413])', pca_sdf_out(session_aud_neurons(pair_i),1000+[sound_onset_ms(1):sound_onset_ms(1)+413])');
        areal_corr_signal_2(count,:) = corr(pca_sdf_out(session_frontal_neurons(pair_i),1000+[sound_onset_ms(2):sound_onset_ms(2)+413])', pca_sdf_out(session_aud_neurons(pair_i),1000+[sound_onset_ms(2):sound_onset_ms(2)+413])');
        areal_corr_signal_3(count,:) = corr(pca_sdf_out(session_frontal_neurons(pair_i),1000+[sound_onset_ms(3):sound_onset_ms(3)+413])', pca_sdf_out(session_aud_neurons(pair_i),1000+[sound_onset_ms(3):sound_onset_ms(3)+413])');
        areal_corr_signal_4(count,:) = corr(pca_sdf_out(session_frontal_neurons(pair_i),1000+[sound_onset_ms(4):sound_onset_ms(4)+413])', pca_sdf_out(session_aud_neurons(pair_i),1000+[sound_onset_ms(4):sound_onset_ms(4)+413])');
        areal_corr_signal_5(count,:) = corr(pca_sdf_out(session_frontal_neurons(pair_i),1000+[sound_onset_ms(5):sound_onset_ms(5)+413])', pca_sdf_out(session_aud_neurons(pair_i),1000+[sound_onset_ms(5):sound_onset_ms(5)+413])');
    end

end

figuren;
histogram(areal_corr_baseline, -1:0.1:1, 'LineStyle', 'none')
histogram(areal_corr_signal_1, -1:0.1:1, 'LineStyle', 'none')
histogram(areal_corr_signal_2, -1:0.1:1, 'LineStyle', 'none')
histogram(areal_corr_signal_3, -1:0.1:1, 'LineStyle', 'none')
histogram(areal_corr_signal_4, -1:0.1:1, 'LineStyle', 'none')
histogram(areal_corr_signal_5, -1:0.1:1, 'LineStyle', 'none')

vline(median(areal_corr_baseline));
vline(median(areal_corr_signal_1));
vline(median(areal_corr_signal_2));
vline(median(areal_corr_signal_3));
vline(median(areal_corr_signal_4));
vline(median(areal_corr_signal_5));