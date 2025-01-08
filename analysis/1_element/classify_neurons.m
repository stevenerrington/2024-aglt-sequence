function neuron_class = classify_neurons(spike_log, element_neuron_table)
    % Define default brain areas if not provided
    brain_areas.auditory = {'R', 'A1', 'RM', 'dSTS'};
    brain_areas.frontal = {'44', '45', 'FOP'};

    % Define indices for neurons in specific brain regions
    auditory_neuron_idx = find(ismember(spike_log.area, brain_areas.auditory));
    frontal_neuron_idx = find(ismember(spike_log.area, brain_areas.frontal));

    % Identify neurons modulated by ROC analysis
    modulated_neurons = find(element_neuron_table.sig_glm == 1);

    % Identify non-modulated neurons
    nonmodulated_neurons = setdiff(1:size(spike_log, 1), modulated_neurons);

    % Find modulated and non-modulated neurons in auditory and frontal areas
    aud_mod_neurons = intersect(modulated_neurons, auditory_neuron_idx);
    frontal_mod_neurons = intersect(modulated_neurons, frontal_neuron_idx);

    aud_nonmod_neurons = intersect(nonmodulated_neurons, auditory_neuron_idx);
    frontal_nonmod_neurons = intersect(nonmodulated_neurons, frontal_neuron_idx);

    % Identify auditory neurons that are positively modulated
    aud_mod_neurons_pos = aud_mod_neurons(strcmp(element_neuron_table.glm_dir_label(aud_mod_neurons), 'pos'));
    % Identify auditory neurons that are negatively modulated
    aud_mod_neurons_neg = aud_mod_neurons(strcmp(element_neuron_table.glm_dir_label(aud_mod_neurons), 'neg'));
    % Identify auditory neurons that are mixed modulated
    aud_mod_neurons_mixed = aud_mod_neurons(strcmp(element_neuron_table.glm_dir_label(aud_mod_neurons), 'mixed'));

    % Identify frontal neurons that are positively modulated
    frontal_mod_neurons_pos = frontal_mod_neurons(strcmp(element_neuron_table.glm_dir_label(frontal_mod_neurons), 'pos'));
    % Identify frontal neurons that are negatively modulated
    frontal_mod_neurons_neg = frontal_mod_neurons(strcmp(element_neuron_table.glm_dir_label(frontal_mod_neurons), 'neg'));
    % Identify frontal neurons that are mixed modulated
    frontal_mod_neurons_mixed = frontal_mod_neurons(strcmp(element_neuron_table.glm_dir_label(frontal_mod_neurons), 'mixed'));

    % Assign results to the output struct
    neuron_class.auditory.all = auditory_neuron_idx;
    neuron_class.auditory.modulated = aud_mod_neurons;
    neuron_class.auditory.nonmodulated = aud_nonmod_neurons;
    neuron_class.auditory.glm_pos = aud_mod_neurons_pos;
    neuron_class.auditory.glm_neg = aud_mod_neurons_neg;
    neuron_class.auditory.glm_mix = aud_mod_neurons_mixed;

    neuron_class.frontal.all = frontal_neuron_idx;
    neuron_class.frontal.modulated = frontal_mod_neurons;
    neuron_class.frontal.nonmodulated = frontal_nonmod_neurons;
    neuron_class.frontal.glm_pos = frontal_mod_neurons_pos;
    neuron_class.frontal.glm_neg = frontal_mod_neurons_neg;
    neuron_class.frontal.glm_mix = frontal_mod_neurons_mixed;
end