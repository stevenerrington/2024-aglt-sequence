neuron_i = 1498;

fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

sdf_in = load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data\spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
event_table_in = load(fullfile(dirs.mat_data,[spike_log.session{neuron_i} '.mat']),'event_table');

sdf_detrended = []; acorr = [];

for trial_i = 1:size(event_table_in.event_table,1)
    sdf_trial = []; sdf_trial = smooth(sdf_in.sdf.sequenceOnset(trial_i,1000+[0:3000]),100)';

    % Linear regression to detrend
    p = polyfit(0:3000, sdf_trial, 1);
    trend = polyval(p, 0:3000);
    sdf_detrended(trial_i, :) = sdf_trial - trend;

    [acorr(trial_i,:), lags] = xcorr(sdf_detrended(trial_i, :), 'coeff');
    sdf_out(trial_i,:) = smooth(sdf_in.sdf.sequenceOnset(trial_i,:),100)';
end


pos_acorr = acorr(:,lags>=0);
mean_acorr = nanmean(pos_acorr);

%% 

trials_in = find(strcmp(event_table_in.event_table.cond_label,'nonviol'));

clear single_unit_fig
% Plot 1: Example neuron raster plot
single_unit_fig(1, 1) = gramm('x', sdf_in.raster.sequenceOnset(trials_in));
single_unit_fig(1, 1).geom_raster('geom', {'point'});
single_unit_fig(1, 1).axe_property('XLim', [-250 3000]);
single_unit_fig(1, 1).axe_property('XTick', [], 'XColor', [1 1 1]);

% Plot 2: Example neuron SDF plot
single_unit_fig(2, 1) = gramm('x', ops.timewin, 'y', sdf_out(trials_in,:));
single_unit_fig(2, 1).stat_summary();
single_unit_fig(2, 1).axe_property('XLim', [-250 3000], 'YLim', [0 30]);

single_unit_fig(1, 1).geom_vline('XIntercept',sound_times);
single_unit_fig(2, 1).geom_vline('XIntercept',sound_times);

single_unit_fig(1, 1).set_layout_options('Position', [0.15 0.75 0.75 0.2], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
single_unit_fig(2, 1).set_layout_options('Position', [0.15 0.2 0.75 0.5], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);

% Create the figure window with specific size and renderer settings
figure('Renderer', 'painters', 'Position', [100 100 500 400]);
single_unit_fig.draw();  % Draw all plots in the figure

%% 
% Plot 2: Example neuron SDF plot
autocorr_example_figure(1, 1) = gramm('x', lags, 'y', acorr(trials_in,:));
autocorr_example_figure(1, 1).geom_line('alpha',0.05);
autocorr_example_figure(1, 1).stat_summary();

autocorr_example_figure(1, 1).axe_property('XLim', [-1400 1400], 'YLim', [-0.5 1]);

% Create the figure window with specific size and renderer settings
figure('Renderer', 'painters', 'Position', [100 100 500 400]);
autocorr_example_figure.draw();  % Draw all plots in the figure
