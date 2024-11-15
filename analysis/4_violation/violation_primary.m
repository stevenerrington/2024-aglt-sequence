clear glm_output encoding_flag encoding_beta z_sdf viol_sdf nonviol_sdf

% Run Violation GLM
% ----------------------------------------------------------------
% Parallel loop to perform GLM analysis for each neuron
parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    cond_inc = [];
    cond_inc = [1 5 14 2 6 15 3 7 13 4 8 16];

    [glm_output{neuron_i}, encoding_flag(neuron_i,:), encoding_beta(neuron_i,:), z_sdf{neuron_i}] = violation_detect_glm(sdf_in, event_table_in)

    viol_sdf(neuron_i,:) = nanmean(z_sdf{neuron_i}(ismember(event_table_in.event_table.cond_value,cond_inc ) & strcmp(event_table_in.event_table.cond_label, 'viol'), :));
    nonviol_sdf(neuron_i,:) = nanmean(z_sdf{neuron_i}(ismember(event_table_in.event_table.cond_value,cond_inc ) & strcmp(event_table_in.event_table.cond_label, 'nonviol'), :));
    
end


% Identify neurons with significant modulation based on GLM results
% ----------------------------------------------------------------
glm_sig_units = find(encoding_flag > 0); % Find neurons with significant encoding in any sound category

pos_glm_sig_units = glm_sig_units(encoding_beta(glm_sig_units) > 0 )
neg_glm_sig_units = glm_sig_units(encoding_beta(glm_sig_units) < 0 )

frontal_viol_neurons = intersect(glm_sig_units,neuron_class.frontal.all);
frontal_viol_neurons_pos = intersect(pos_glm_sig_units,neuron_class.frontal.all);
frontal_viol_neurons_neg = intersect(neg_glm_sig_units,neuron_class.frontal.all);

auditory_viol_neurons = intersect(glm_sig_units,neuron_class.auditory.all);
auditory_viol_neurons_pos = intersect(pos_glm_sig_units,neuron_class.auditory.all);
auditory_viol_neurons_neg = intersect(neg_glm_sig_units,neuron_class.auditory.all);


% Plot figure
% ----------------------------------------------------------------
dev_violation_glm_example


% Decoding analysis (population)
% ----------------------------------------------------------------

clear violation_class_data nonviolation_class_data
violation_class_data = viol_sdf(frontal_viol_neurons_neg,:);
nonviolation_class_data = nonviol_sdf(frontal_viol_neurons_neg,:);

%% Accuracy x time

obs_data = {violation_class_data, nonviolation_class_data};
obs_labels = {'viol','nonviol'};

accuracy = run_multicat_categorization(obs_data, obs_labels);

