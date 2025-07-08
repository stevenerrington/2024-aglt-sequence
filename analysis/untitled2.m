% Clear previous PCA and nonviolent SDF data
clear pca_* nonviol_sdf

% Loop through each neuron using parallel for-loop for efficiency
for neuron_i = 1:size(spike_log,1)
    % Display current progress
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load event table for the corresponding session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviolent' trials
    nonviol_sdf = [];
    nonviol_sdf = sdf_in.sdf.sequenceOnset( ...
        strcmp(event_table_in.event_table.cond_label, 'nonviol') & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms), :);

    % Calculate mean and standard deviation of baseline firing rate (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, 800:1000)));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, 800:1000)));

    % Compute normalized and smoothed SDF for all nonviolent trials
    pca_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);

    % Generate a shuffled version of the SDF (temporal shuffling within trials)
    pca_sdf_out_shuffled(neuron_i,:) = smooth((nanmean(nonviol_sdf(:,randperm(size(nonviol_sdf, 2)))) - baseline_fr_mean) ./ baseline_fr_std, 50);

    % Compute normalized and smoothed SDFs for each sequence condition
    pca_sdf_out_seq1(neuron_i,:) = smooth( ...
        (nanmean(sdf_in.sdf.sequenceOnset( ...
        (event_table_in.event_table.cond_value == 1 | event_table_in.event_table.cond_value == 5) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms), :)) - baseline_fr_mean) ./ baseline_fr_std, 50);

    pca_sdf_out_seq2(neuron_i,:) = smooth( ...
        (nanmean(sdf_in.sdf.sequenceOnset( ...
        (event_table_in.event_table.cond_value == 2 | event_table_in.event_table.cond_value == 6) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms), :)) - baseline_fr_mean) ./ baseline_fr_std, 50);

    pca_sdf_out_seq3(neuron_i,:) = smooth( ...
        (nanmean(sdf_in.sdf.sequenceOnset( ...
        (event_table_in.event_table.cond_value == 3 | event_table_in.event_table.cond_value == 7) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms), :)) - baseline_fr_mean) ./ baseline_fr_std, 50);

    pca_sdf_out_seq4(neuron_i,:) = smooth( ...
        (nanmean(sdf_in.sdf.sequenceOnset( ...
        (event_table_in.event_table.cond_value == 4 | event_table_in.event_table.cond_value == 8) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms), :)) - baseline_fr_mean) ./ baseline_fr_std, 50);

    n_trials_seq(neuron_i,:) = [sum((event_table_in.event_table.cond_value == 1 | event_table_in.event_table.cond_value == 5) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms)),...
        sum((event_table_in.event_table.cond_value == 2 | event_table_in.event_table.cond_value == 6) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms)),...
        sum((event_table_in.event_table.cond_value == 3 | event_table_in.event_table.cond_value == 7) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms)),...
        sum((event_table_in.event_table.cond_value == 4 | event_table_in.event_table.cond_value == 8) & ...
        ~isnan(event_table_in.event_table.rewardOnset_ms))];
end

% Clean trials with less than 10 trials
nonvalid_seq_neurons = find(any(n_trials_seq < 10,2));
for non_valid_i = 1:length(nonvalid_seq_neurons)
    pca_sdf_out(nonvalid_seq_neurons(non_valid_i),:) = nan(1,6001);
    pca_sdf_out_seq1(nonvalid_seq_neurons(non_valid_i),:) = nan(1,6001);
    pca_sdf_out_seq2(nonvalid_seq_neurons(non_valid_i),:) = nan(1,6001);
    pca_sdf_out_seq3(nonvalid_seq_neurons(non_valid_i),:) = nan(1,6001);
    pca_sdf_out_seq4(nonvalid_seq_neurons(non_valid_i),:) = nan(1,6001);
end


%% Run PCA
pc_out_auditory = perform_pca_and_plot(auditory_neuron_idx, pca_sdf_out);
pc_out_frontal = perform_pca_and_plot(frontal_neuron_idx, pca_sdf_out);

%% Plot PCA
% Time indices corresponding to sound events (in ms)
sound_times = [0, 563, 1126, 1689, 2252];
sound_times_idx = sound_times + 100; % Adjust index for plotting

% Plot 3D PCA trajectory for auditory neurons
figuren('Renderer', 'painters', 'Position', [200 200 1000 400]); hold on;
subplot(1, 2, 1); hold on
color_line3(pc_out_auditory.obs.pcs(:,1), pc_out_auditory.obs.pcs(:,2), pc_out_auditory.obs.pcs(:,3), -100:2750, 'LineWidth', 2);
view(326.4106,22.1996);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); grid on
scatter3(pc_out_auditory.obs.pcs(sound_times_idx,1), pc_out_auditory.obs.pcs(sound_times_idx,2), pc_out_auditory.obs.pcs(sound_times_idx,3), 100, [0 0 0], '^', 'filled');

% Plot 3D PCA trajectory for frontal neurons
subplot(1, 2, 2); hold on
color_line3(pc_out_frontal.obs.pcs(:,1), pc_out_frontal.obs.pcs(:,2), pc_out_frontal.obs.pcs(:,3), -100:2750, 'LineWidth', 2);
view(326.4106,22.1996);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); grid on
scatter3(pc_out_frontal.obs.pcs(sound_times_idx,1), pc_out_frontal.obs.pcs(sound_times_idx,2), pc_out_frontal.obs.pcs(sound_times_idx,3), 100, [0 0 0], '^', 'filled');


%% Starting point ]]] IN PROGRESS [[[
% Time indices corresponding to sound events (in ms)
sound_times = [0, 563, 1126, 1689, 2252];
sound_times_idx = sound_times + 100; % Adjust index for plotting
clear onset_distances_*
for boot_i = 1:nboot
    fprintf('Iteration: %i of %i \n', boot_i, nboot);

    clear pc_out_auditory pc_out_frontal
    pc_out_auditory = perform_pca_and_plot(randsample(auditory_neuron_idx,100), pca_sdf_out);
    pc_out_frontal = perform_pca_and_plot(randsample(frontal_neuron_idx,100), pca_sdf_out);

    clear pca_traj1 pca_traj2 pca1_z pca2_z
    pca_traj1 = [pc_out_auditory.obs.pcs(:,1), pc_out_auditory.obs.pcs(:,2), pc_out_auditory.obs.pcs(:,3)];
    pca_traj2 = [pc_out_frontal.obs.pcs(:,1), pc_out_frontal.obs.pcs(:,2), pc_out_frontal.obs.pcs(:,3)];

    pca1_z = zscore(pca_traj1);
    pca2_z = zscore(pca_traj2);

    for element_i = 3:5
        onset_distances_auditory(boot_i,element_i-2) = euclidean_distance_3d(pca1_z(sound_times_idx(2),:), pca1_z(sound_times_idx(element_i),:));
        onset_distances_frontal(boot_i,element_i-2) = euclidean_distance_3d(pca2_z(sound_times_idx(2),:), pca2_z(sound_times_idx(element_i),:));
    end

end

plot_distance_data = [nanmean(onset_distances_auditory,2); nanmean(onset_distances_frontal,2)];
plot_distance_label = [repmat({'Auditory'}, length(onset_distances_auditory),1); repmat({'Frontal'}, length(onset_distances_frontal),1)];

% Plot bootstrap classification accuracy
figure('Renderer', 'painters', 'Position', [100 100 500 500]);
clear plot_bootstrap_distance
plot_bootstrap_distance(1,1) = gramm('x', plot_distance_label, 'y', plot_distance_data, 'color', plot_distance_label);
plot_bootstrap_distance(1,1).stat_summary('geom', {'bar', 'black_errorbar'}, 'dodge', 1, 'width', 0.5);
plot_bootstrap_distance(1,1).axe_property('YLim', [2 3]);
plot_bootstrap_distance.draw;



plot_distanceSep_data_auditory = [onset_distances_auditory(:,1); onset_distances_auditory(:,2); onset_distances_auditory(:,3); onset_distances_auditory(:,4)];
plot_distanceSep_data_frontal = [onset_distances_frontal(:,1); onset_distances_frontal(:,2); onset_distances_frontal(:,3); onset_distances_frontal(:,4)];
plot_distanceSep_label = [repmat({'1 -> 2'},1000,1); repmat({'1 -> 3'},1000,1); repmat({'1 -> 4'},1000,1); repmat({'1 -> 5'},1000,1)];
plot_distanceSep_area = [repmat({'Auditory'},4000,1); repmat({'Frontal'},4000,1)];

% Plot bootstrap classification accuracy
figure('Renderer', 'painters', 'Position', [100 100 500 500]);
clear plot_bootstrap_distanceSep
plot_bootstrap_distanceSep(1,1) = gramm('x', [plot_distanceSep_label; plot_distanceSep_label], 'y', [plot_distanceSep_data_auditory; plot_distanceSep_data_frontal], 'color', plot_distanceSep_area);
plot_bootstrap_distanceSep(1,1).stat_summary('geom', {'point', 'line', 'black_errorbar'}, 'width', 0.5);
plot_bootstrap_distanceSep(1,1).axe_property('YLim', [1 5]);
plot_bootstrap_distanceSep.draw;



%% Trajectory comparison
