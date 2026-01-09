%% ------------------------------------------------------------------------
% Select significant neurons (auditory + frontal)
% -------------------------------------------------------------------------
sig_neurons = [
    neuron_class.auditory.all;
    neuron_class.frontal.all
];

%% ------------------------------------------------------------------------
% Prepare SDF data
% -------------------------------------------------------------------------
sdf_raw = normFR_in.norm_fr_soundAll(sig_neurons,:);   % [nNeurons x nTime]

% Preallocate
sdf_cluster = zeros(size(sdf_raw));
sdf_plot    = zeros(size(sdf_raw));

% Smooth SDFs
for i = 1:size(sdf_raw,1)
    sdf_cluster(i,:) = smooth(sdf_raw(i,:), 10);   % light smoothing (clustering)
    sdf_plot(i,:)    = smooth(sdf_raw(i,:), 50);   % heavy smoothing (plotting)
end

%% ------------------------------------------------------------------------
% Run MDS + clustering pipeline
% -------------------------------------------------------------------------
mds_neuron_pipeline(sdf_cluster)

% (Optional) recompute nonmetric MDS stress explicitly
[Y, stress] = mdscale(mds_results.D, 5);
fprintf('Nonmetric MDS stress = %.3f\n', stress);

%% ------------------------------------------------------------------------
% Plot cluster-average motifs
% -------------------------------------------------------------------------
cluster_order = [3 9 1 6 13 7 5 10 11 8 4 12 0 2];

ylim_clusters = {
    [-2 8], [-2 8], [-2 8], [-2 8], [-2 6], [-2 4], [-2 6],...
    [-8 2], [-6 2], [-8 2], [-4 6], [-6 2], [-1 1], [-6 2]
};

figuren('Renderer','painters','Position',[100 100 1400 350]);

for ii = 1:numel(cluster_order)
    k = cluster_order(ii);

    if k == 0
        subplot(2,7,ii); hold on;
        continue
    else
    subplot(2,7,ii); hold on;

    plot(-200:800, mds_results.cluster_mean(k,:), 'LineWidth',1.5);

    vline([0 563], 'k-');    % stimulus window
    vline(413, 'k--');       % event of interest
    axis square
    xlim([-150 600]);
    ylim(ylim_clusters{ii});
    title(sprintf('Cluster %d', k));
    end
end

%% ------------------------------------------------------------------------
% Count neurons per cluster
% -------------------------------------------------------------------------

for cluster_i = 1:13
    n_per_cluster(cluster_i,1) = cluster_i;
    n_per_cluster(cluster_i,2) = sum(mds_results.cluster_idx == cluster_i);
end

n_per_cluster(cluster_order(cluster_order > 0),2)

%% ------------------------------------------------------------------------
% Group clusters by response type
% -------------------------------------------------------------------------
facilitated_clusters  = [3 9 1 6 13 7];
suppressed_clusters   = [10 11 8 4 12];
facilitated_ramping   = 5;
suppressed_ramping    = 2;

cluster_groups = {
    facilitated_clusters
    suppressed_clusters
    facilitated_ramping
    suppressed_ramping
};

group_labels = {
    'Facilitated'
    'Suppressed'
    'Facilitated ramping'
    'Suppressed ramping'
};

% Map neurons to grouped cluster types
cluster_idx_grouped = cell(numel(cluster_groups),1);

for g = 1:numel(cluster_groups)
    cluster_idx_grouped{g} = sig_neurons( ...
        ismember(mds_results.cluster_idx, cluster_groups{g}) );
end

%% ------------------------------------------------------------------------
% Count auditory vs frontal neurons per group
% -------------------------------------------------------------------------
n_neurons_cluster = zeros(numel(cluster_groups),2);

for g = 1:numel(cluster_groups)
    n_neurons_cluster(g,1) = numel(intersect( ...
        cluster_idx_grouped{g}, auditory_neuron_idx));

    n_neurons_cluster(g,2) = numel(intersect( ...
        cluster_idx_grouped{g}, frontal_neuron_idx));

    cluster_neurons.auditory{g} = intersect( ...
        cluster_idx_grouped{g}, auditory_neuron_idx);
    
    cluster_neurons.frontal{g} = intersect( ...
        cluster_idx_grouped{g}, frontal_neuron_idx);

end

% Normalize within each area
n_neurons_cluster(:,1) = n_neurons_cluster(:,1) ./ sum(n_neurons_cluster(:,1));
n_neurons_cluster(:,2) = n_neurons_cluster(:,2) ./ sum(n_neurons_cluster(:,2));

%% ------------------------------------------------------------------------
% Plot stacked bar chart (area comparison)
% -------------------------------------------------------------------------
figure('Renderer','painters','Position',[200 200 500 350]); hold on;

bar(n_neurons_cluster','stacked','LineStyle','none');

xticks([1 2]);
xticklabels({'Auditory','Frontal'});
ylabel('Proportion of neurons');

legend(group_labels,'Location','best');
box off;


% Rows = cortex (auditory, frontal)
% Columns = neuron category (fac, sup, fac ramp, sup ramp)
% Run the chi-square test
[chi2_stat, df, p_value, expected] = chi2_independence(n_neurons_cluster');

fprintf('Chi-square = %.2f, df = %d, p = %.4f\n', chi2_stat, df, p_value);
disp('Expected counts under independence:');
disp(expected);


%%
% Example transient, suppression, and ramping neurons

clu_exp_neurons = [1359 1280 721 1352];
plot_ylim_example = {[-2 10], [-2 14], [30 70], [30 50]};
plot_startxpos = [0.1 0.35 0.6 0.85]-0.05;
order = [2 1 3 4];

% Reorder cluster-specific variables
clu_exp_neurons = clu_exp_neurons(order);
plot_ylim_example = plot_ylim_example(order);

% Load color palette
colorscale = cbrewer('qual','Dark2',6);
color_pal = struct();
for clu_i = 1:4
    color_pal.(['clu' int2str(clu_i)]) = colorscale(clu_i,:);
end

clear example_* single_unit_fig

% Loop through clusters to plot raster, example SDF, and population SDF
for neuron_i = 1:4
    example_neuron_idx = clu_exp_neurons(neuron_i);
    xlim_vals = [-180 600];

    % Select valid trials (non-violative, non-baseline)
    example_valid_idx = strcmp(string(sdf_soundAlign_data{example_neuron_idx}(:,5)), 'nonviol') & ...
                        ~strcmp(string(sdf_soundAlign_data{example_neuron_idx}(:,3)), 'Baseline');

    % Extract SDF and raster
    example_sdf_in = cell2mat(sdf_soundAlign_data{example_neuron_idx}(example_valid_idx,1));
    example_raster_in = sdf_soundAlign_data{example_neuron_idx}(example_valid_idx,2);

    % Smooth individual trial SDF
    for trial_i = 1:size(example_sdf_in,1)
        example_sdf_in(trial_i,:) = smooth(example_sdf_in(trial_i,:),50)';
    end

    % Raster plot
    single_unit_fig(1, neuron_i) = gramm('x', example_raster_in);
    single_unit_fig(1, neuron_i).geom_raster('geom', {'point'});
    single_unit_fig(1, neuron_i).axe_property('XLim', xlim_vals);
    single_unit_fig(1, neuron_i).set_color_options('map', color_pal.(['clu' int2str(neuron_i)]));

    % Example neuron SDF
    single_unit_fig(2, neuron_i) = gramm('x', ops.sound_sdf_window, 'y', example_sdf_in);
    single_unit_fig(2, neuron_i).stat_summary();
    single_unit_fig(2, neuron_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_example{neuron_i});
    single_unit_fig(2, neuron_i).set_color_options('map', color_pal.(['clu' int2str(neuron_i)]));

    % Layout for subplots
    single_unit_fig(1, neuron_i).set_layout_options('Position', [plot_startxpos(neuron_i) 0.7 0.18 0.2], ...
        'legend', false, 'margin_height', [0 0], 'margin_width', [0 0], 'redraw', false);
    single_unit_fig(2, neuron_i).set_layout_options('Position', [plot_startxpos(neuron_i) 0.1 0.18 0.5], ...
        'margin_height', [0 0], 'margin_width', [0 0], 'redraw', false);

    % Customize plot appearance
    single_unit_fig(1, neuron_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(2, neuron_i).geom_vline('xintercept', [0 413 563], 'style', {'k-', 'k--', 'k-'});
end

figure('Renderer', 'painters', 'Position', [100 100 1200 350]);
single_unit_fig.draw();




