clear seq_sdf_out
% Loop through each neuron
parfor neuron_i = 1:size(spike_log,1)
    % Display progress for the current neuron
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviolent' condition
    nonviol_sdf = [];
    nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label, 'nonviol') & ~isnan(event_table_in.event_table.rewardOnset_ms), :);
    nonviol_rwd_sdf = sdf_in.sdf.reward(strcmp(event_table_in.event_table.cond_label, 'nonviol') & ~isnan(event_table_in.event_table.rewardOnset_ms), :);

    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 50);
    seq_sdf_out_rwd(neuron_i,:) = smooth((nanmean(nonviol_rwd_sdf) - baseline_fr_mean) ./ baseline_fr_std, 50);
end


clusterid = 'clu2';

figuren('Renderer', 'painters', 'Position', [100 100 1800 300]); hold on;
a1 = subplot(1,4,[1 2 3]); hold on
plot(-1000:5000, nanmean(seq_sdf_out(intersect(auditory_neuron_idx,neuron_class.cluster_idx.(clusterid)),:)))
plot(-1000:5000, nanmean(seq_sdf_out(intersect(frontal_neuron_idx,neuron_class.cluster_idx.(clusterid)),:)))
xlim([-200 3000]); box off
vline([0 563 1126 1689 2252],'k-')
vline([413 976 1539 2102 2665],'k--')

a2 = subplot(1,4,[4]); hold on
plot(-1000:5000, nanmean(seq_sdf_out_rwd(intersect(auditory_neuron_idx,neuron_class.cluster_idx.(clusterid)),:)))
plot(-1000:5000, nanmean(seq_sdf_out_rwd(intersect(frontal_neuron_idx,neuron_class.cluster_idx.(clusterid)),:)))
xlim([-800 800]); box off

linkaxes([a1 a2],'y')


%%

figuren('Renderer', 'painters', 'Position', [100 100 500 300]); hold on;
a1 = subplot(1,2,1); hold on
plot(-1000:5000, nanmean(seq_sdf_out_rwd(intersect(auditory_neuron_idx, intersect(neuron_class.nhp.troy,neuron_class.cluster_idx.(clusterid))),:)))
plot(-1000:5000, nanmean(seq_sdf_out_rwd(intersect(auditory_neuron_idx, intersect(neuron_class.nhp.walt,neuron_class.cluster_idx.(clusterid))),:)))
xlim([-800 800]); box off

a1 = subplot(1,2,2); hold on
plot(-1000:5000, nanmean(seq_sdf_out_rwd(intersect(frontal_neuron_idx, intersect(neuron_class.nhp.troy,neuron_class.cluster_idx.(clusterid))),:)))
plot(-1000:5000, nanmean(seq_sdf_out_rwd(intersect(frontal_neuron_idx, intersect(neuron_class.nhp.walt,neuron_class.cluster_idx.(clusterid))),:)))
xlim([-800 800]); box off




%% Example neuron
