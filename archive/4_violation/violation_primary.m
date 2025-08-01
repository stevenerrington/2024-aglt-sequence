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

    %[glm_output{neuron_i}, encoding_flag(neuron_i,:), encoding_beta(neuron_i,:), z_sdf{neuron_i}] = violation_detect_glm(sdf_in, event_table_in)


    [roc_neuron_out_ROC, roc_neuron_out_p, sdf_viol_out, sdf_nonviol_out] = violation_detect_ROC(spike_log, dirs)

    viol_trials = []; viol_trials = ismember(event_table_in.event_table.cond_value,cond_inc ) & strcmp(event_table_in.event_table.cond_label, 'viol') & ~isnan(event_table_in.event_table.rewardOnset_ms);
    nonviol_trials = []; nonviol_trials = ismember(event_table_in.event_table.cond_value,cond_inc ) & strcmp(event_table_in.event_table.cond_label, 'nonviol') & ~isnan(event_table_in.event_table.rewardOnset_ms);

    viol_sdf(neuron_i,:) = nanmean(z_sdf{neuron_i}(viol_trials, :));
    nonviol_sdf(neuron_i,:) = nanmean(z_sdf{neuron_i}(nonviol_trials, :));
    
end


% Identify neurons with significant modulation based on GLM results
% ----------------------------------------------------------------
glm_sig_units = find(encoding_flag > 0); % Find neurons with significant encoding in any sound category

pos_glm_sig_units = glm_sig_units(encoding_beta(glm_sig_units) > 0 );
neg_glm_sig_units = glm_sig_units(encoding_beta(glm_sig_units) < 0 );

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


%%



%% Accuracy x time

obs_data = {violation_class_data, nonviolation_class_data};
obs_labels = {'viol','nonviol'};

accuracy = run_multicat_categorization(obs_data, obs_labels);

%%
%% Distance/similarity analysis
% Define the time window and extract indices for the time window from 'ops'
timewin = [-50:1:500];  % Set the time window from -100 ms to 2750 ms
timewin_idx = find(ismember([-200:800], timewin));  % Get the indices of this time window from the 'ops' structure

% Store the spike density function (SDF) signals for each sequence in the specified time window
signal_in_frontal = {};  % Initialize a cell array to store the signal data
signal_in_frontal = {viol_sdf(neuron_class.frontal.all, timewin_idx),  % SDF for sequence 1 in the current area
    nonviol_sdf(neuron_class.frontal.all, timewin_idx)}; % SDF for sequence 4

% Store the spike density function (SDF) signals for each sequence in the specified time window
signal_in_auditory = {};  % Initialize a cell array to store the signal data
signal_in_auditory = {viol_sdf(neuron_class.auditory.all, timewin_idx),  % SDF for sequence 1 in the current area
    nonviol_sdf(neuron_class.auditory.all, timewin_idx)}; % SDF for sequence 4


% Perform cross-condition PCA analysis to extract principal components
[pc_out_frontal, pc_shuf_out_frontal] = get_xcond_pca(signal_in_frontal);  % Perform PCA on the input signals
[pc_out_auditory, pc_shuf_out_auditory] = get_xcond_pca(signal_in_auditory);  % Perform PCA on the input signals


pc1a = pc_out_auditory{1}(1,:);
pc2a = pc_out_auditory{1}(2,:);
pc3a = pc_out_auditory{1}(3,:);

pc1b = pc_out_auditory{2}(1,:);
pc2b = pc_out_auditory{2}(2,:);
pc3b = pc_out_auditory{2}(3,:);


figuren('Renderer', 'painters', 'Position', [100 207 634 593]);
subplot(3,1,1); hold on; box off
plot(timewin, pc1a);plot(timewin, pc1b)

subplot(3,1,2); hold on; box off
plot(timewin, pc2a);plot(timewin, pc2b)

subplot(3,1,3); hold on; box off
plot(timewin, pc3a);plot(timewin, pc3b)

legend({'Violation','Nonviolation'})

pc1a = pc_out_frontal{1}(1,:);
pc2a = pc_out_frontal{1}(2,:);
pc3a = pc_out_frontal{1}(3,:);

pc1b = pc_out_frontal{2}(1,:);
pc2b = pc_out_frontal{2}(2,:);
pc3b = pc_out_frontal{2}(3,:);


figuren('Renderer', 'painters', 'Position', [100 207 634 593]);
subplot(3,1,1); hold on; box off
plot(timewin, pc1a);plot(timewin, pc1b)

subplot(3,1,2); hold on; box off
plot(timewin, pc2a);plot(timewin, pc2b)

subplot(3,1,3); hold on; box off
plot(timewin, pc3a);plot(timewin, pc3b)

legend({'Violation','Nonviolation'})