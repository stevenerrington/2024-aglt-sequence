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
% Plot frotnal neuron PCA trajectory
% ------------------------------------------------------------
figure('Renderer', 'painters', 'Position', [200 200 1000 400]); hold on;

subplot(1,1,1); hold on

plot(smooth(pc_out_frontal.obs.pcs(:,2), smooth_factor), ...
    smooth(pc_out_frontal.obs.pcs(:,3), smooth_factor), 'LineWidth', 1.5);
xlabel('PC2'); ylabel('PC3'); grid on; axis square;

scatter(pc_out_frontal.obs.pcs(sound_times_idx,2), ...
         pc_out_frontal.obs.pcs(sound_times_idx,3), ...
         100, [0 0 0], 'o', 'filled');
title('Frontal Neurons');


%%
x = smooth(pc_out_frontal.obs.pcs(:,2), smooth_factor);
y = smooth(pc_out_frontal.obs.pcs(:,3), smooth_factor);

x_2 = smooth(pc_out_auditory.obs.pcs(:,2), smooth_factor);
y_2 = smooth(pc_out_auditory.obs.pcs(:,3), smooth_factor);

t = pc_out_frontal.window; % time vector for coloring

figure('Renderer', 'painters', 'Position', [200 200 1000 400]); hold on;
subplot(1,2,1); hold on
scatter(x, y, 20, t, 'filled');
colormap(parula);
colorbar;
xlabel('PC2'); ylabel('PC3'); grid on; axis square;
title('Frontal Neurons');
scatter(pc_out_frontal.obs.pcs(sound_times_idx,2), ...
        pc_out_frontal.obs.pcs(sound_times_idx,3), ...
        200, t(sound_times_idx), ...
        'filled');

subplot(1,2,2); hold on
scatter(x_2, y_2, 20, t, 'filled');
colormap(parula);
colorbar;
xlabel('PC2'); ylabel('PC3'); grid on; axis square;
title('Auditory Neurons');
scatter(pc_out_auditory.obs.pcs(sound_times_idx,2), ...
        pc_out_auditory.obs.pcs(sound_times_idx,3), ...
        200, t(sound_times_idx), ...
        'filled');


%%

x_frontal = smooth(pc_out_frontal.obs.pcs(:,1), smooth_factor);
x_auditory = smooth(pc_out_auditory.obs.pcs(:,1), smooth_factor);
t = pc_out_frontal.window;

figure('Renderer', 'painters', 'Position', [200 200 600 400]); hold on;
plot(t, x_frontal./max(abs(x_frontal)))
plot(t, x_auditory./max(abs(x_auditory)))
ylim([-1.2 1.2]); grid on
xlim([min(t) max(t)])