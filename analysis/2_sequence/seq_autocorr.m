rng(1,"twister"); % Set random number generator seed for reproducibility using 'twister' algorithm
acorr_time = 1000+[0:2665]; % Define time window for autocorrelation (ms)

%%
sound_times = [0, 563, 1126, 1689, 2252]; % Define sound onset times (ms) relative to some reference
acorr_neuron = [];
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

    seq_sdf_out_temp = seq_sdf_out(neuron_i,:);
    [acorr_neuron(neuron_i,:), ~] = xcorr(seq_sdf_out_temp(:,acorr_time), 'coeff');

end

%% Run cluster average autocorrelation

% Compute population-averaged autocorrelation for auditory and frontal neurons
[acorr_aud_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.auditory.all,acorr_time)), 'coeff');
[acorr_frontal_pos, lags] = xcorr(nanmean(seq_sdf_out(neuron_class.frontal.all,acorr_time)), 'coeff');

% Compute autocorrelation for each cluster (1 to 6) within auditory and frontal neuron classes
for clu_i = 1:6
    [acorr_auditory(clu_i,:), lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all),acorr_time)), 'coeff');
    [acorr_frontal(clu_i,:), lags] = xcorr(nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.frontal.all),acorr_time)), 'coeff');
end

%% Get peak estimates 

for clu_i = 1:6
    clear peak_idx_aud peak_idx_frontal
    [~,peak_idx_aud] = findpeaks(acorr_auditory(clu_i,:),'MinPeakProminence',0.4);
    [~,peak_idx_frontal] = findpeaks(acorr_frontal(clu_i,:),'MinPeakProminence',0.4);

    try
        peak_clu(clu_i,1) = lags(peak_idx_aud(find(lags(peak_idx_aud) > 0,1)));
    catch
        peak_clu(clu_i,1) = NaN;
    end

    try
        peak_clu(clu_i,2) = lags(peak_idx_frontal(find(lags(peak_idx_frontal) > 0,1)));
    catch
        peak_clu(clu_i,2) = NaN;
    end
end


%% Figure: plot seq sdf
figure_xlim = [-250 2665]; % X-axis limits for plots
figure_linewidth = 1; % Line width for plots
figuren('Renderer', 'painters', 'Position', [21,700,1371,229]); hold on

% Y-axis limits for each cluster (auditory and frontal)
plot_ylim_clu_i = {[-0.5 2.5], [-0.25 0.5];...
    [-0.2 0.8], [-0.2 0.2];...
    [-0.4 0.6], [-0.2 0.2];...
    [-0.6 0.4], [-0.2 0.4];...
    [-0.6 0.2], [-0.4 0.2];...
    [-0.2 0.6], [-0.2 0.6];};

% Plot average SDFs for each cluster (auditory top row, frontal bottom row)
for clu_i = 1:6
    subplot(2,6,clu_i); hold on
    plot(-1000:5000, nanmean(seq_sdf_out(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all),:)),'LineWidth',figure_linewidth,'Color',color_pal.(['clu' int2str(clu_i)]))
    xlim(figure_xlim); ylim(plot_ylim_clu_i{clu_i,1})
    vline([sound_onset_ms], 'k-') % Vertical line at sound onset
    vline([sound_onset_ms]+413, 'k--') % Vertical line offset by 413 ms

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
xlim([-1500 1500]); vline([-1126 -563 0 563 1126],'k--') % Mark sound time intervals
title('Auditory')

subplot(1,2,2); hold on; box off
for clu_i = 1:6
    plot(lags,acorr_frontal(clu_i,:),'color',color_pal.(['clu' int2str(clu_i)]),'LineWidth',1);
end
xlim([-1500 1500]); vline([-1126 -563 0 563 1126],'k--')
title('Frontal')

%%
example_seq_neurons = [neuron_class.cluster_idx.clu6(30), neuron_class.cluster_idx.clu1(6)]; % Select example neurons
acorr_time = 1000+[0:2665]; % Time window for autocorrelation

% Process selected example neurons
for neuron_i = 1:length(example_seq_neurons)
    example_sequence_neuron = example_seq_neurons(neuron_i);

    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{example_sequence_neuron} '_' spike_log.unitDSP{example_sequence_neuron} '.mat']));

    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{example_sequence_neuron} '.mat']), 'event_table');

    nonviol_trials = find(strcmp(event_table_in.event_table.cond_label, 'nonviol')); % Identify nonviol trials

    % Extract spike density function (SDF) for 'nonviol' condition
    nonviol_sdf = []; nonviol_raster = [];
    nonviol_sdf = sdf_in.sdf.sequenceOnset(nonviol_trials, :);
    nonviol_raster = sdf_in.raster.sequenceOnset(nonviol_trials);

    % Smooth each trial's SDF with a 100-point smoothing window
    for trial_i = 1:size(nonviol_sdf,1)
        nonviol_sdf(trial_i,:) = smooth(nonviol_sdf(trial_i,:),100)';
    end

    % Compute autocorrelation of the mean SDF
    acorr = []; lags = [];
    [acorr, lags] = xcorr(nanmean(nonviol_sdf(:,acorr_time)), 'coeff');

    % Store outputs for example neuron
    seq_sdfEx_out{neuron_i} = nonviol_sdf;
    seq_raster_out{neuron_i} = nonviol_raster;
    seq_acorr_out{neuron_i} = acorr;
    seq_lags_out{neuron_i} = lags;
end

% Get example neuron index for the current cluster
xlim_vals = [-500 2800];  % Set X-axis limits for plots
ylim_vals = {[20 70], [0 20]}; % Set Y-axis limits for example SDF plots

clear seq_example_fig
seq_example_fig_xstart = [0.1 0.6]; % X start positions for subplots

% Plot raster and SDF for each example neuron
for neuron_i = 1:length(example_seq_neurons)

    % Plot 1: Raster plot for example neuron
    seq_example_fig(1, neuron_i) = gramm('x', seq_raster_out{neuron_i});
    seq_example_fig(1, neuron_i).geom_raster('geom', {'point'});
    seq_example_fig(1, neuron_i).axe_property('XLim', xlim_vals);

    seq_example_fig(1,neuron_i).set_layout_options...
    ('Position',[seq_example_fig_xstart(neuron_i) 0.8 0.35 0.15],... % Set subplot position
    'legend',false,...
    'margin_height',[0.00 0.00],... % Set margins
    'margin_width',[0.0 0.00],...
    'redraw',false);

    % Plot 2: SDF plot for example neuron (mean ± SEM)
    seq_example_fig(2, neuron_i) = gramm('x', -1000:5000, 'y', seq_sdfEx_out{neuron_i});
    seq_example_fig(2, neuron_i).stat_summary();
    seq_example_fig(2, neuron_i).axe_property('XLim', xlim_vals, 'YLim', ylim_vals{neuron_i});
    seq_example_fig(2, neuron_i).geom_vline('xintercept',[0 516 1126 1689 2252],'style','k-');

    seq_example_fig(2,neuron_i).set_layout_options...
    ('Position',[seq_example_fig_xstart(neuron_i) 0.2 0.35 0.5],... % Set subplot position
    'legend',false,...
    'margin_height',[0.00 0.00],...
    'margin_width',[0.0 0.00],...
    'redraw',false);
end

% Draw figure with raster and SDF plots for example neurons
figure('Renderer', 'painters', 'Position', [117,500,916,373]);
seq_example_fig.draw;

% Plot autocorrelation for the example neurons
figuren('Renderer', 'painters', 'Position', [1041,572,400,300]);
plot(seq_lags_out{1}, seq_acorr_out{1},'LineWidth',1.5)
plot(seq_lags_out{2}, seq_acorr_out{2},'LineWidth',1.5)
vline([sound_onset_ms],'k--')
vline([-sound_onset_ms],'k--')
axis square; xlim([-1000 1000]); ylim([0 1])
