rng(1,"twister");

%%
sound_times = [0, 563, 1126, 1689, 2252];

% Loop through each neuron
parfor neuron_i = 1:size(spike_log,1)
    % Display progress for the current neuron
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviol' condition
    nonviol_sdf = [];
    nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label, 'nonviol'), :);

    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    seq_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 100);
end


%% Run autocorrelation

acorr_time = 1000+[0:2665];

[acorr_aud_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.all,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.all,acorr_time)), 'coeff');

for clu_i = 1:6
    [acorr_auditory(clu_i,:), lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all),acorr_time)), 'coeff');
    [acorr_frontal(clu_i,:), lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.frontal.all),acorr_time)), 'coeff');

end


%% Figure: plot seq sdf
figure_xlim = [-250 2665];
figure_linewidth = 1;
figuren('Renderer', 'painters', 'Position', [21,700,1371,229]); hold on

plot_ylim_clu_i = {[-0.5 2.5], [-0.25 0.5];...
    [-0.2 0.8], [-0.2 0.2];...
    [-0.4 0.6], [-0.2 0.2];...
    [-0.6 0.4], [-0.2 0.4];...
    [-0.6 0.2], [-0.4 0.2];...
    [-0.2 0.6], [-0.2 0.6];};


for clu_i = 1:6
    subplot(2,6,clu_i); hold on
    plot(-1000:5000, nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all),:)),'LineWidth',figure_linewidth,'Color',color_pal.(['clu' int2str(clu_i)]))
    xlim(figure_xlim); ylim(plot_ylim_clu_i{clu_i,1})
    vline([sound_onset_ms], 'k-')
    vline([sound_onset_ms]+413, 'k--')

    subplot(2,6,clu_i+6); hold on
    plot(-1000:5000, nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.frontal.all),:)),'LineWidth',figure_linewidth,'Color',color_pal.(['clu' int2str(clu_i)]))
    xlim(figure_xlim); ylim(plot_ylim_clu_i{clu_i,2})
    vline([sound_onset_ms], 'k-')
    vline([sound_onset_ms]+413, 'k--')
end

%%

figuren('Renderer', 'painters', 'Position', [1317,696,476,226]); hold on
subplot(1,2,1); hold on; box off
for clu_i = 1:6
    plot(lags,acorr_auditory(clu_i,:),'color',color_pal.(['clu' int2str(clu_i)]),'LineWidth',1);
end
xlim([-1500 1500]); vline([-1126 -563 0 563 1126],'k--')
title('Auditory')

subplot(1,2,2); hold on; box off
for clu_i = 1:6
    plot(lags,acorr_frontal(clu_i,:),'color',color_pal.(['clu' int2str(clu_i)]),'LineWidth',1);
end
xlim([-1500 1500]); vline([-1126 -563 0 563 1126],'k--')
title('Frontal')

