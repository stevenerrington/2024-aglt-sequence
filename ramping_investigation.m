

figure;
subplot(2,2,1)
plot(ops.timewin, nanmean(pca_sdf_out(cluster_neurons.frontal{3},:)))
vline(sound_onset_ms); vline(sound_onset_ms+413,'r-'); xlim([-250 3500])
title('Frontal - pos')

subplot(2,2,3)
plot(ops.timewin, nanmean(pca_sdf_out(cluster_neurons.frontal{4},:)))
vline(sound_onset_ms,'r-'); vline(sound_onset_ms+413,'r--'); xlim([-250 3500])
title('Frontal - neg')

subplot(2,2,2)
plot(ops.timewin, nanmean(pca_sdf_out(cluster_neurons.auditory{3},:)))
vline(sound_onset_ms,'r-'); vline(sound_onset_ms+413,'r--'); xlim([-250 3500])
title('Auditory - pos')

subplot(2,2,4)
plot(ops.timewin, nanmean(pca_sdf_out(cluster_neurons.auditory{4},:)))
vline(sound_onset_ms,'r-'); vline(sound_onset_ms+413,'r--'); xlim([-250 3500])
title('Auditory - neg')


%%

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

    nonviol_sdf = sdf_in.sdf.reward(nonviol_mask, :);

    % ------------------------------------------------------------
    % Baseline normalization (–200 ms to 0 ms relative window)
    % ------------------------------------------------------------
    baseline_window = 4500:5000; % Adjust if sampling window differs
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, baseline_window)));
    baseline_fr_std  = nanstd(nanmean(nonviol_sdf(:, baseline_window)));

    % ------------------------------------------------------------
    % Compute smoothed & normalized SDF (mean over nonviolent trials)
    % ------------------------------------------------------------
    pca_sdf_out_rwd(neuron_i,:) = smooth( ...
        (nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);


end

%%

figure;
subplot(2,2,1)
plot(ops.timewin, nanmean(pca_sdf_out_rwd(cluster_neurons.frontal{3},:)))
xlim([-1000 1000])
title('Frontal - pos')

subplot(2,2,3)
plot(ops.timewin, nanmean(pca_sdf_out_rwd(cluster_neurons.frontal{4},:)))
xlim([-1000 1000])
title('Frontal - neg')

subplot(2,2,2)
plot(ops.timewin, nanmean(pca_sdf_out_rwd(cluster_neurons.auditory{3},:)))
xlim([-1000 1000])
title('Auditory - pos')

subplot(2,2,4)
plot(ops.timewin, nanmean(pca_sdf_out_rwd(cluster_neurons.auditory{4},:)))
xlim([-1000 1000])
title('Auditory - neg')

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
        isnan(event_table_in.event_table.rewardOnset_ms) &...
        ~isnan(event_table_in.event_table.stimulusOnset_ms);

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
    pca_sdf_out_nonrwd(neuron_i,:) = smooth( ...
        (nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);

end

for neuron_i = 1:size(spike_log,1)
    if any(isinf(pca_sdf_out_nonrwd(neuron_i,:)))
        pca_sdf_out_nonrwd(neuron_i,:) = nan(1,6001);
    end
end

%%


figure;
subplot(2,2,1)
plot(ops.timewin, nanmean(pca_sdf_out_nonrwd(cluster_neurons.frontal{3},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Frontal - pos')

subplot(2,2,3)
plot(ops.timewin, nanmean(pca_sdf_out_nonrwd(cluster_neurons.frontal{4},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Frontal - neg')

subplot(2,2,2)
plot(ops.timewin, nanmean(pca_sdf_out_nonrwd(cluster_neurons.auditory{3},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Auditory - pos')

subplot(2,2,4)
plot(ops.timewin, nanmean(pca_sdf_out_nonrwd(cluster_neurons.auditory{4},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Auditory - neg')

%%
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
    nonviol_mask = strcmp(event_table_in.event_table.cond_label, 'viol') & ...
        ~isnan(event_table_in.event_table.stimulusOnset_ms);

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
    pca_sdf_out_viol(neuron_i,:) = smooth( ...
        (nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);

end

for neuron_i = 1:size(spike_log,1)
    if any(isinf(pca_sdf_out_nonrwd(neuron_i,:)))
        pca_sdf_out_viol(neuron_i,:) = nan(1,6001);
    end
end


%%

figure;
subplot(2,2,1)
plot(ops.timewin, nanmean(pca_sdf_out_viol(cluster_neurons.frontal{3},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Frontal - pos')

subplot(2,2,3)
plot(ops.timewin, nanmean(pca_sdf_out_viol(cluster_neurons.frontal{4},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Frontal - neg')

subplot(2,2,2)
plot(ops.timewin, nanmean(pca_sdf_out_viol(cluster_neurons.auditory{3},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Auditory - pos')

subplot(2,2,4)
plot(ops.timewin, nanmean(pca_sdf_out_viol(cluster_neurons.auditory{4},:)))
vline(sound_onset_ms); xlim([-250 3500])
title('Auditory - neg')

%%
session_i = 10;
datafile = ephysLog.session{session_i};
data_in = load(fullfile(dirs.mat_data,datafile))
fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), datafile)
