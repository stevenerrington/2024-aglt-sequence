% Define indices of example neurons for each cluster
clu_exp_neurons = [18 9 8 12 2 30];

% Set y-axis limits for the example neuron plots (SDF)
plot_ylim_example = {[2 8], [0 10], [0 50], [0 15], [6 14], [30 50]};

% Set y-axis limits for the population SDF plots
plot_ylim_pop = {[-4 4], [-5 20], [-10 5], [-2 8], [-10 4], [-3 6]};

% Define horizontal positioning for each cluster plot in the figure layout
plot_startxpos = [0.1 0.25 0.4 0.55 0.7 0.85];

% Define the desired order of clusters for plotting
order = [2 1 4 3 5 6];

% Reorder cluster-specific variables based on the specified order
clu_exp_neurons = clu_exp_neurons(order);
plot_ylim_example = plot_ylim_example(order);
plot_ylim_pop = plot_ylim_pop(order);

% Load color palette for cluster visualizations
colorscale = cbrewer('qual','Dark2',6);

% Initialize color palette structure for each cluster
color_pal = struct();
for clu_i = 1:6
    color_pal.(['clu' int2str(clu_i)]) = colorscale(clu_i,:);
end

% Clear previously stored variables starting with 'example_'
clear example_* single_unit_fig

% Loop over each of the 6 clusters to extract and plot data
for cluster_i = 1:6
    % Get example neuron index for the current cluster
    example_neuron_idx = neuron_class.cluster_idx.(['clu' int2str(cluster_i)])(clu_exp_neurons(cluster_i));
    xlim_vals = [-180 600];  % Set X-axis limits for plots

    % Identify trials that are non-violative and not baseline
    example_valid_idx = strcmp(string(sdf_soundAlign_data{example_neuron_idx}(:,5)) ,'nonviol') & ~strcmp(string(sdf_soundAlign_data{example_neuron_idx}(:,3)) ,'Baseline');

    % Extract spike density functions (SDF) for valid trials
    example_sdf_in = cell2mat(sdf_soundAlign_data{example_neuron_idx}(example_valid_idx, 1));
    
    % Extract raster data for valid trials
    example_raster_in = sdf_soundAlign_data{example_neuron_idx}(example_valid_idx, 2);

    % Smooth each trial's SDF
    for trial_i = 1:size(example_sdf_in,1)
        example_sdf_in(trial_i,:) = smooth(example_sdf_in(trial_i,:),100)';
    end

    % Plot 1: Raster plot for example neuron
    single_unit_fig(1, cluster_i) = gramm('x', example_raster_in);
    single_unit_fig(1, cluster_i).geom_raster('geom', {'point'});
    single_unit_fig(1, cluster_i).axe_property('XLim', xlim_vals);
    single_unit_fig(1, cluster_i).set_color_options('map',color_pal.(['clu' int2str(cluster_i)]));

    % Plot 2: SDF plot for example neuron (mean ± SEM)
    single_unit_fig(2, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', example_sdf_in);
    single_unit_fig(2, cluster_i).stat_summary();
    single_unit_fig(2, cluster_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_example{cluster_i});
    single_unit_fig(2, cluster_i).set_color_options('map',color_pal.(['clu' int2str(cluster_i)]));

    % Plot 3: Population SDF for current cluster
    single_unit_fig(3, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', inputSDF_plot{1}(neuron_class.(['cluster' int2str(cluster_i)]),:));
    single_unit_fig(3, cluster_i).stat_summary();
    single_unit_fig(3, cluster_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_pop{cluster_i});
    single_unit_fig(3, cluster_i).set_color_options('map',color_pal.(['clu' int2str(cluster_i)]));

    % Define subplot layout for raster, example SDF, and population SDF
    single_unit_fig(1, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.8 0.1 0.1], 'legend', false, 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(2, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.5 0.1 0.25], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(3, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.2 0.1 0.25], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);

    % Customize plot appearance (remove ticks from upper plots, add event lines)
    single_unit_fig(1, cluster_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(2, cluster_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(3, cluster_i).geom_vline('xintercept', 0, 'style', 'k-');
    single_unit_fig(3, cluster_i).geom_vline('xintercept', 413, 'style', 'k--');
    single_unit_fig(3, cluster_i).geom_vline('xintercept', 563, 'style', 'k-');
end

% Create figure and draw all gramm plots
figure('Renderer', 'painters', 'Position', [100 100 1200 500]);
single_unit_fig.draw();

%% P(neurons) per area

% Define cluster labels for plotting purposes
cluster_labels = {'onset','facilitated','biphasic', 'offset', 'suppressed', 'ramping'};

% Calculate total neuron count per brain area
sum(n_area)

% Normalize neuron counts to get proportions for each area
p_n_area(:,1) = n_area(:,1)./sum(n_area(:,1));  % Auditory
p_n_area(:,2) = n_area(:,2)./sum(n_area(:,2));  % Frontal

% Create bar plot showing proportion of each cluster per area
figuren('Renderer', 'painters', 'Position', [200 200 500 350]); hold on;
bh = bar(p_n_area','stacked','LineStyle','none');  % Transpose to plot area as x-axis
xticks([1 2]); xticklabels({'Auditory','Frontal'})  % Label x-axis
legend(cluster_labels,'Location','eastoutside')  % Add legend for cluster types

% Apply cluster-specific colors to the bars
set(bh, 'FaceColor', 'Flat')
for clu_i = 1:6
    bh(clu_i).CData = color_pal.(['clu' int2str(clu_i)]);
end
