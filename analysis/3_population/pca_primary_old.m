% Clear previous PCA and nonviolent SDF data
clear pca_* nonviol_sdf

% Initialize an array to keep track of neurons that caused errors
error_neuron = [];

% Loop through each neuron using parallel for-loop for efficiency
parfor neuron_i = 1:size(spike_log,1)
    try
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

    catch
        % Log any neuron that caused an error
        error_neuron = [error_neuron; neuron_i];
    end
end

% Remove the neuron with the maximum value at time index 4250 (possible outlier)
pca_sdf_out(find(pca_sdf_out(:,4250) == max(pca_sdf_out(:,4250)),1),:) = nan(1,6001);
pca_sdf_out_seq1(find(pca_sdf_out_seq1(:,4250) == max(pca_sdf_out_seq1(:,4250)),1),:) = nan(1,6001);
pca_sdf_out_seq2(find(pca_sdf_out_seq2(:,4250) == max(pca_sdf_out_seq2(:,4250)),1),:) = nan(1,6001);
pca_sdf_out_seq3(find(pca_sdf_out_seq3(:,4250) == max(pca_sdf_out_seq3(:,4250)),1),:) = nan(1,6001);
pca_sdf_out_seq4(find(pca_sdf_out_seq4(:,4250) == max(pca_sdf_out_seq4(:,4250)),1),:) = nan(1,6001);

%% Perform PCA and Plot for Each Brain Region

% Perform PCA for auditory neurons
pc_out_auditory = perform_pca_and_plot(auditory_neuron_idx, pca_sdf_out);

% Perform PCA for frontal neurons
pc_out_frontal = perform_pca_and_plot(frontal_neuron_idx, pca_sdf_out);

% Plot cumulative explained variance of PCs
figuren('Renderer', 'painters', 'Position', [200 200 400 400]); hold on
plot(cumsum(pc_out_auditory.obs.var_exp),'Color', color_pal.auditory_clu(4,:),'LineWidth',1.5)
plot(cumsum(pc_out_frontal.obs.var_exp),'Color', color_pal.frontal_clu(4,:),'LineWidth',1.5)
xlim([0 50]); xlabel('Principal Component'); ylabel('')
hline(80,'k') % Horizontal line at 80% variance explained

% Display how much variance is explained by the first 6 PCs
cumsum(pc_out_auditory.obs.var_exp(1:6));
cumsum(pc_out_frontal.obs.var_exp(1:6));

%% 3D PCA Trajectory Visualization

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

%% Calculate distance between trajectory start points

% Number of bootstrap iterations
nboot = 100;

% Number of neurons to sample in each bootstrap
nbootsamples = 20;

% Initialize arrays to store bootstrapped neuron indices
boot_neurons = [];

% Generate bootstrap samples for auditory and frontal neurons
for boot_i = 1:nboot
    % Randomly sample 'nbootsamples' auditory neurons (with replacement) and sort them
    aud_boot_neurons(:, boot_i) = sort(auditory_neuron_idx(randperm(length(auditory_neuron_idx), nbootsamples)));
    
    % Randomly sample 'nbootsamples' frontal neurons (with replacement) and sort them
    frontal_boot_neurons(:, boot_i) = sort(frontal_neuron_idx(randperm(length(frontal_neuron_idx), nbootsamples)));
end

% Loop over each bootstrap iteration
for boot_i = 1:nboot
    % Clear previous PCA output variables
    clear pc_out*

    % Perform PCA on the current bootstrap sample of auditory neurons
    pc_out_auditory = perform_pca_and_plot(aud_boot_neurons(:, boot_i), pca_sdf_out);

    % Perform PCA on the current bootstrap sample of frontal neurons
    pc_out_frontal = perform_pca_and_plot(frontal_boot_neurons(:, boot_i), pca_sdf_out);

    % Initialize arrays to hold PCA coordinates at specific time points (element onset)
    auditory_elementonset_xyz = [];
    frontal_elementonset_xyz = [];

    % Extract PC1–PC3 values at sound time indices for auditory neurons
    auditory_elementonset_xyz = [ ...
        pc_out_auditory.obs.pcs(sound_times_idx, 1), ...
        pc_out_auditory.obs.pcs(sound_times_idx, 2), ...
        pc_out_auditory.obs.pcs(sound_times_idx, 3)];

    % Extract PC1–PC3 values at sound time indices for frontal neurons
    frontal_elementonset_xyz = [ ...
        pc_out_frontal.obs.pcs(sound_times_idx, 1), ...
        pc_out_frontal.obs.pcs(sound_times_idx, 2), ...
        pc_out_frontal.obs.pcs(sound_times_idx, 3)];

    % Compute pairwise Euclidean distances between each pair of element onset points
    for ele_i = 1:5
        for ele_j = 1:5
            % Distance between elements i and j in auditory PCA space
            distance_auditory(ele_i, ele_j, boot_i) = euclidean_distance_3d( ...
                auditory_elementonset_xyz(ele_i,:), auditory_elementonset_xyz(ele_j,:));
            
            % Distance between elements i and j in frontal PCA space
            distance_frontal(ele_i, ele_j, boot_i) = euclidean_distance_3d( ...
                frontal_elementonset_xyz(ele_i,:), frontal_elementonset_xyz(ele_j,:));
        end
    end
end
