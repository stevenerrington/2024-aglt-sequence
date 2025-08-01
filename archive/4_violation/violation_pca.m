rng(1,"twister");
 
% Loop through each neuron
parfor neuron_i = 1:size(spike_log,1)
    % Display progress for the current neuron
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviolation' condition
    nonviol_sdf = [];
    nonviol_sdf = sdf_in.sdf.violation(strcmp(event_table_in.event_table.cond_label, 'nonviol') & ~isnan(event_table_in.event_table.rewardOnset_ms), :);
    
    % Extract spike density function (SDF) for 'violation' condition
    viol_sdf = [];
    viol_sdf = sdf_in.sdf.violation(strcmp(event_table_in.event_table.cond_label, 'viol') & ~isnan(event_table_in.event_table.rewardOnset_ms), :);

    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean_nonviol = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std_nonviol = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_mean_viol = nanmean(nanmean(viol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std_viol = nanstd(nanmean(viol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    seq_nonviol_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean_nonviol) ./ baseline_fr_std_nonviol, 100);
    seq_viol_out(neuron_i,:) = smooth((nanmean(viol_sdf) - baseline_fr_mean_viol) ./ baseline_fr_std_viol, 100);

    if nanmean(seq_nonviol_out(neuron_i,:)) > 50
        seq_nonviol_out(neuron_i,:) = nan(1,length(ops.timewin))
    end
    if nanmean(seq_viol_out(neuron_i,:)) > 50
        seq_viol_out(neuron_i,:) = nan(1,length(ops.timewin))
    end    

end



%%
figuren('Renderer', 'painters', 'Position', [100 100 1500 500]); 
subplot(1,2,1); hold on
plot(ops.timewin, nanmean(seq_nonviol_out(frontal_neuron_idx,:)))
plot(ops.timewin, nanmean(seq_viol_out(frontal_neuron_idx,:)))
xlim([-563 2252]); vline(sound_onset_ms,'k-')

subplot(1,2,2); hold on
plot(ops.timewin, nanmean(seq_nonviol_out(auditory_neuron_idx,:)))
plot(ops.timewin, nanmean(seq_viol_out(auditory_neuron_idx,:)))
xlim([-563 2252]); vline(sound_onset_ms,'k-')


%%
timewin = [-100:5:600];  % Set the time window from -100 ms to 2750 ms
timewin_idx = find(ismember(ops.timewin, timewin));  % Get the indices of this time window from the 'ops' structure

clear signal_in* pc_out*
% Store the spike density function (SDF) signals for each sequence in the specified time window
signal_in_frontal = {};  % Initialize a cell array to store the signal data
signal_in_frontal = {seq_nonviol_out(frontal_neuron_idx,timewin_idx),...
    seq_viol_out(frontal_neuron_idx,timewin_idx)}; % SDF for sequence 4

signal_in_auditory = {};  % Initialize a cell array to store the signal data
signal_in_auditory = {seq_nonviol_out(auditory_neuron_idx,timewin_idx),...
    seq_viol_out(auditory_neuron_idx,timewin_idx)}; % SDF for sequence 4

% Perform cross-condition PCA analysis to extract principal components
[pc_out_frontal, pc_shuf_out_frontal] = get_xcond_pca(signal_in_frontal);  % Perform PCA on the input signals
[pc_out_auditory, pc_shuf_out_auditory] = get_xcond_pca(signal_in_auditory);  % Perform PCA on the input signals

figuren('Renderer', 'painters', 'Position', [100 100 1500 500]); 
subplot(1,2,1); hold on
zero_idx = find(timewin == 0);
plot3(pc_out_frontal{1}(1,:), pc_out_frontal{1}(2,:), pc_out_frontal{1}(3,:))
plot3(pc_out_frontal{2}(1,:), pc_out_frontal{2}(2,:), pc_out_frontal{2}(3,:))
scatter3(pc_out_frontal{1}(1,zero_idx), pc_out_frontal{1}(2,zero_idx), pc_out_frontal{1}(3,zero_idx),'k','filled')
scatter3(pc_out_frontal{2}(1,zero_idx), pc_out_frontal{2}(2,zero_idx), pc_out_frontal{2}(3,zero_idx),'k','filled')

subplot(1,2,2); hold on
zero_idx = find(timewin == 0);
plot3(signal_in_auditory{1}(1,:), signal_in_auditory{1}(3,:), signal_in_auditory{1}(4,:))
plot3(signal_in_auditory{2}(1,:), signal_in_auditory{2}(2,:), signal_in_auditory{2}(3,:))
scatter3(signal_in_auditory{1}(1,zero_idx), signal_in_auditory{1}(2,zero_idx), signal_in_auditory{1}(3,zero_idx),'k','filled')
scatter3(signal_in_auditory{2}(1,zero_idx), signal_in_auditory{2}(2,zero_idx), signal_in_auditory{2}(3,zero_idx),'k','filled')