%% ============================================================
%  PCA Analysis of Neuronal Spike Density Functions (SDF)
%  ------------------------------------------------------------
%  This script processes single-unit spike data, computes
%  normalized SDFs, performs PCA, and visualizes trajectories
%  for auditory and frontal neuron populations.
% ============================================================

%% Setup
% Clear previous PCA variables and specific SDF data
clear pca_* nonviol_sdf

%% ------------------------------------------------------------
% Loop through each neuron (using parfor if available)
% ------------------------------------------------------------
for neuron_i = 1:size(spike_log,1)

    % Display progress in console
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load spike data and event table for this neuron
    sdf_in = load(fullfile(dirs.root, 'data', 'spike', ...
        [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
    event_table_in = load(fullfile(dirs.mat_data, ...
        [spike_log.session{neuron_i} '.mat']), 'event_table');

    % ------------------------------------------------------------
    % Extract nonviolation SDFs (only rewarded trials)
    % ------------------------------------------------------------
    nonviol_mask = strcmp(event_table_in.event_table.cond_label, 'nonviol') & ...
                   ~isnan(event_table_in.event_table.rewardOnset_ms);

    nonviol_sdf = sdf_in.sdf.sequenceOnset(nonviol_mask, :);

    % ------------------------------------------------------------
    % Baseline normalization (–200 ms to 0 ms relative window)
    % ------------------------------------------------------------
    baseline_window = 800:1000; % Adjust if sampling window differs
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, baseline_window)));
    baseline_fr_std  = nanstd(nanmean(nonviol_sdf(:, baseline_window)));

    % ------------------------------------------------------------
    % Compute smoothed & normalized SDF (mean over nonviolent trials)
    % ------------------------------------------------------------
    pca_sdf_out(neuron_i,:) = smooth( ...
        (nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);

    % Shuffled control (temporal shuffling within trials)
    pca_sdf_out_shuffled(neuron_i,:) = smooth( ...
        (nanmean(nonviol_sdf(:, randperm(size(nonviol_sdf, 2)))) - baseline_fr_mean) ./ baseline_fr_std, 50);

    % ------------------------------------------------------------
    % Condition-specific SDFs by sequence number (1–4, 5–8)
    % ------------------------------------------------------------
    for seq = 1:4
        seq_mask = (event_table_in.event_table.cond_value == seq | ...
                    event_table_in.event_table.cond_value == seq + 4) & ...
                    ~isnan(event_table_in.event_table.rewardOnset_ms);

        seq_sdf = sdf_in.sdf.sequenceOnset(seq_mask, :);
        eval(sprintf('pca_sdf_out_seq%i(neuron_i,:) = smooth((nanmean(seq_sdf) - baseline_fr_mean) ./ baseline_fr_std, 50);', seq))

        n_trials_seq(neuron_i, seq) = sum(seq_mask);
    end

end

%% ------------------------------------------------------------
% Remove neurons with insufficient trials (< 10 per sequence)
% ------------------------------------------------------------
nonvalid_seq_neurons = find(any(n_trials_seq < 10, 2));

for i = 1:length(nonvalid_seq_neurons)
    idx = nonvalid_seq_neurons(i);
    pca_sdf_out(idx,:) = nan(1,6001);
    for seq = 1:4
        eval(sprintf('pca_sdf_out_seq%i(idx,:) = nan(1,6001);', seq));
    end
end


%% ============================================================
% Run PCA on auditory and frontal neurons
% ============================================================
pc_out_auditory = perform_pca_and_plot(auditory_neuron_idx, pca_sdf_out);
pc_out_frontal  = perform_pca_and_plot(frontal_neuron_idx,  pca_sdf_out);

% Display cumulative variance explained (first 6 PCs)
disp(cumsum(pc_out_auditory.obs.var_exp(1:6)))
disp(cumsum(pc_out_frontal.obs.var_exp(1:6)))


%% ============================================================
% Plot 3D PCA Trajectories
% ============================================================

% Define sound onset times (ms)
sound_times = [0, 563, 1126, 1689, 2252];
closest_idx  = zeros(size(sound_times));

for i = 1:length(sound_times)
    [~, closest_idx(i)] = min(abs(pc_out_auditory.window - sound_times(i)));
end

sound_times_idx = closest_idx;
smooth_factor = 1;

% Define color maps
aud_color = plasma; 
frontal_color = plasma;
aud_color = max(min(aud_color, 1), 0);
frontal_color = max(min(frontal_color, 1), 0);

% ------------------------------------------------------------
% Plot auditory neuron PCA trajectory
% ------------------------------------------------------------
figure('Renderer', 'painters', 'Position', [200 200 1000 400]); hold on;

subplot(1,2,1); hold on
color_line3( ...
    smooth(pc_out_auditory.obs.pcs(:,1), smooth_factor), ...
    smooth(pc_out_auditory.obs.pcs(:,2), smooth_factor), ...
    smooth(pc_out_auditory.obs.pcs(:,3), smooth_factor), ...
    pc_out_auditory.window, 'LineWidth', 1.5, 'colormap', aud_color);
view(401, 25);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); grid on; axis square;
scatter3(pc_out_auditory.obs.pcs(sound_times_idx,1), ...
         pc_out_auditory.obs.pcs(sound_times_idx,2), ...
         pc_out_auditory.obs.pcs(sound_times_idx,3), ...
         100, [0 0 0], 'o', 'filled');
title('Auditory Neurons');

% ------------------------------------------------------------
% Plot frontal neuron PCA trajectory
% ------------------------------------------------------------
subplot(1,2,2); hold on
color_line3( ...
    smooth(pc_out_frontal.obs.pcs(:,1), smooth_factor), ...
    smooth(pc_out_frontal.obs.pcs(:,2), smooth_factor), ...
    smooth(pc_out_frontal.obs.pcs(:,3), smooth_factor), ...
    pc_out_frontal.window, 'LineWidth', 1.5, 'colormap', frontal_color);
view(401, 25);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); grid on; axis square;
scatter3(pc_out_frontal.obs.pcs(sound_times_idx,1), ...
         pc_out_frontal.obs.pcs(sound_times_idx,2), ...
         pc_out_frontal.obs.pcs(sound_times_idx,3), ...
         100, [0 0 0], 'o', 'filled');
title('Frontal Neurons');


%% ============================================================
% Bootstrap Analysis of Trajectory Metrics
% ============================================================

nboot = 1000;
fprintf('\nBootstrapping %i iterations...\n', nboot);

for boot_i = 1:nboot
    fprintf('Iteration: %i of %i \n', boot_i, nboot);

    % Bootstrap resampling of neurons (with replacement)
    pc_out_auditory = perform_pca_and_plot(randsample(auditory_neuron_idx, 100, true), pca_sdf_out);
    pc_out_frontal  = perform_pca_and_plot(randsample(frontal_neuron_idx,  100, true), pca_sdf_out);

    % Extract normalized PCA trajectories
    pca1_z = zscore(pc_out_auditory.obs.pcs(:,1:3), 1, 'all');
    pca2_z = zscore(pc_out_frontal.obs.pcs(:,1:3),  1, 'all');

    % Quantify trajectory dynamics
    metrics_auditory = quantify_trajectory_dynamics(pca1_z, -100:2750);
    metrics_frontal  = quantify_trajectory_dynamics(pca2_z, -100:2750);

    % Store results
    traj_pathlength(boot_i,:)       = [metrics_auditory.path_length,       metrics_frontal.path_length];
    traj_curvature(boot_i,:)        = [nanmean(metrics_auditory.curvature), nanmean(metrics_frontal.curvature)];
    traj_net_displacement(boot_i,:) = [metrics_auditory.net_displacement,  metrics_frontal.net_displacement];
    traj_efficiency(boot_i,:)       = [metrics_auditory.efficiency,        metrics_frontal.efficiency];
    traj_R_g(boot_i,:)              = [metrics_auditory.R_g,               metrics_frontal.R_g];
    traj_tortuosity(boot_i,:)       = [metrics_auditory.tortuosity,        metrics_frontal.tortuosity];
    traj_cosinesimilarity(boot_i,:) = [metrics_auditory.mean_cosine_similarity, metrics_frontal.mean_cosine_similarity];

    
    clear distance
    for seq_i = 1:4
        distance(1,seq_i) = norm(mean(pca1_z(sound_times_idx(seq_i) + [-2:1:2],:)) - mean(pca1_z(sound_times_idx(seq_i+1) + [-2:1:2],:)));
        distance(2,seq_i) = norm(mean(pca2_z(sound_times_idx(seq_i) + [-2:1:2],:)) - mean(pca2_z(sound_times_idx(seq_i+1) + [-2:1:2],:)));
    end
    
    traj_intraelementDist(boot_i,:) = [mean(distance(1,:)) mean(distance(2,:))];

end


%% ------------------------------------------------------------
% Compare element onset distance (Auditory vs. Frontal)
% ------------------------------------------------------------




%% ------------------------------------------------------------
% Compare path length (Start–End Distance)
% ------------------------------------------------------------
plot_startend_data_mu = [traj_pathlength(:,1); traj_pathlength(:,2)];
plot_startend_label   = [repmat({'Auditory'}, nboot,1); repmat({'Frontal'}, nboot,1)];

figure('Renderer', 'painters', 'Position', [100 100 200 300]);
g = gramm('x', plot_startend_label, 'y', plot_startend_data_mu, 'color', plot_startend_label);
g.geom_swarm('alpha',0.5,'point_size',0.5);
g.stat_summary('dodge',0.7,'geom',{'black_point','black_errorbar'},'type','quartile');
g.axe_property('YLim', [25 175]);
g.set_names('X','Area','Y','Path Length (a.u.)');
g.no_legend;
g.set_title('Start–End Distance');
g.draw;

% Confidence interval and p-value
onset_diff = traj_pathlength(:,2) - traj_pathlength(:,1);
ci = prctile(onset_diff, [2.5 97.5]);
p_val = 2 * min(mean(onset_diff >= 0), mean(onset_diff <= 0));
disp(p_val)










%% ============================================================
% Linearity of Frontal PC1 (Regression Analysis)
% ============================================================
pc_out_frontal = perform_pca_and_plot(frontal_neuron_idx, pca_sdf_out);
mdl = fitlm(pc_out_frontal.window, pc_out_frontal.obs.pcs(:,1));

figure; 
plot(mdl);
title('Linearity of Frontal PC1');
xlabel('Time (ms)');
ylabel('PC1 Value');
grid on;
