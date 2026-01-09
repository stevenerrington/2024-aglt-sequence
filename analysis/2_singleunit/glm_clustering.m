%% ===============================
%% 1. Identify Significant Neurons
%% ===============================

% Combine auditory and frontal neurons into a single list
sig_neurons = [neuron_class.auditory.all; neuron_class.frontal.all];

%% ===============================
%% 2. Prepare Spike Density Function (SDF) Data
%% ===============================

clear inputSDF

% Extract normalized firing rates for all sound trials
inputSDF = {normFR_in.norm_fr_soundAll(sig_neurons,:)};

% Smooth SDF for clustering and plotting
for neuron_i = 1:size(inputSDF{1},1)
    inputSDF{1}(neuron_i,:) = smooth(inputSDF{1}(neuron_i,:),25);   % Light smoothing for clustering
    inputSDF_plot{1}(neuron_i,:) = smooth(inputSDF{1}(neuron_i,:),25); % Heavier smoothing for plotting
end

% Split SDF by sound conditions
inputSDF_split = {
    normFR_in.norm_fr_soundA(sig_neurons,:),
    normFR_in.norm_fr_soundC(sig_neurons,:),
    normFR_in.norm_fr_soundD(sig_neurons,:),
    normFR_in.norm_fr_soundF(sig_neurons,:),
    normFR_in.norm_fr_soundG(sig_neurons,:)
};

% Define time windows for SDF and epoch analysis
sdfTimes = {[-200:800], [-200:800], [-200:800], [-200:800], [-200:800]};
sdfEpoch = {[-50:563], [-50:563], [-50:563], [-50:563], [-50:563]};

% Define color mapping for clustering
colorMapping = [1, 1, 1, 1, 1];

%% ===============================
%% 3. Consensus Clustering
%% ===============================

[sortIDs, idxDist, raw, respSumStruct, rawLink, myK] = ...
    consensusCluster(inputSDF, sdfTimes, '-e', sdfEpoch, '-ei', colorMapping, '-er', sdfEpoch, '-mn', 25, '-c', 0.25);

nClusters_manual = myK;
clusterNeurons = [];

% Assign neuron indices to clusters
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:, nClusters_manual) == cluster_i);
end

%% ===============================
%% 4. Visualize Clustering
%% ===============================

% Plot dendrogram
figuren('Renderer', 'painters', 'Position', [100 100 500 400]);
subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink, 0, 'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h, rawLink, sortIDs(:, nClusters_manual))
set(gca,'YTick',[],'YLim',[1 length(rawLink)]); 
xlabel('Similarity')

% Plot similarity matrix
subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic, ir) = raw(ir, ic); % Mirror lower triangle
    end
end
imagesc(raw(outPerm, outPerm));
CT = cbrewer('seq', 'PuBuGn', 100); CT(CT>1) = 1; CT(CT<0) = 0;
colormap(flipud(CT));
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim', [-1 1])

% Define color palette for cluster SDF plots
plot_lineclust_color = turbo(nClusters_manual);

%% ===============================
%% 5. Cluster Evaluation: Average SDF
%% ===============================

figuren('Renderer', 'painters', 'Position', [100 100 850 400]);

for cluster_i = 1:nClusters_manual
    subplot(2,5,cluster_i); hold on
    plot(sdfTimes{1}, nanmean(inputSDF_plot{1}(clusterNeurons{cluster_i},:),1), ...
        'Color', plot_lineclust_color(cluster_i,:));
    vline([0 563], 'k-');   % Stimulus window
    vline([413], 'k--');    % Event of interest
    xlim([-150 600])
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end

%% ===============================
%% 6. Manual Cluster Curation
%% ===============================

% Combine all clustered neurons
neuron_class.clusterAll = [];
for cluster_i = 1:nClusters_manual
    neuron_class.(['cluster' int2str(cluster_i)]) = clusterNeurons{cluster_i};
    neuron_class.clusterAll = [neuron_class.clusterAll; clusterNeurons{cluster_i}];
end
neuron_class.clusterAll = sort(neuron_class.clusterAll);

% Identify unclustered neurons
neuron_class.unclustered = find(~ismember(1:size(sig_neurons), neuron_class.clusterAll));

%% ===============================
%% 8. Assign Cluster Indices and Compute Area-Specific Counts
%% ===============================

for clu_i = 1:6
    neuron_class.cluster_idx.(['clu' int2str(clu_i)]) = sig_neurons(neuron_class.(['cluster' int2str(clu_i)]));
    neuron_class.cluster_idx.(['clu' int2str(clu_i) '_aud']) = ...
        intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]), auditory_neuron_idx);
    neuron_class.cluster_idx.(['clu' int2str(clu_i) '_frontal']) = ...
        intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]), frontal_neuron_idx);

    n_clu_area(clu_i,1) = length(neuron_class.cluster_idx.(['clu' int2str(clu_i) '_aud']));
    n_clu_area(clu_i,2) = length(neuron_class.cluster_idx.(['clu' int2str(clu_i) '_frontal']));
end

% Calculate proportion of neurons per area
p_clu_area(:,1) = n_clu_area(:,1)./length(neuron_class.auditory.all);
p_clu_area(:,2) = n_clu_area(:,2)./length(neuron_class.frontal.all);

% Identify unclustered neurons
cluster_cells = arrayfun(@(x) neuron_class.cluster_idx.(['clu' num2str(x)]), 1:6, 'UniformOutput', false);
clustered_neurons = sort(vertcat(cluster_cells{:}));
unclustered_neurons = sig_neurons(~ismember(sig_neurons, clustered_neurons));

%% ===============================
%% 9. Example Neuron Plots per Cluster
%% ===============================

clu_exp_neurons = [29 70 8 12 32 30];
plot_ylim_example = {[0 12], [0 10], [0 50], [0 15], [0 5], [30 50]};
plot_ylim_pop = {[-2 10], [-5 5], [-2 5], [-6 4], [-4 2], [-2 4]};
plot_startxpos = [0.1 0.25 0.4 0.55 0.7 0.85];
order = 1:6;

% Reorder cluster-specific variables
clu_exp_neurons = clu_exp_neurons(order);
plot_ylim_example = plot_ylim_example(order);
plot_ylim_pop = plot_ylim_pop(order);

% Load color palette
colorscale = cbrewer('qual','Dark2',6);
color_pal = struct();
for clu_i = 1:6
    color_pal.(['clu' int2str(clu_i)]) = colorscale(clu_i,:);
end

clear example_* single_unit_fig

% Loop through clusters to plot raster, example SDF, and population SDF
for cluster_i = 1:6
    example_neuron_idx = neuron_class.cluster_idx.(['clu' int2str(cluster_i)])(clu_exp_neurons(cluster_i));
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
    single_unit_fig(1, cluster_i) = gramm('x', example_raster_in);
    single_unit_fig(1, cluster_i).geom_raster('geom', {'point'});
    single_unit_fig(1, cluster_i).axe_property('XLim', xlim_vals);
    single_unit_fig(1, cluster_i).set_color_options('map', color_pal.(['clu' int2str(cluster_i)]));

    % Example neuron SDF
    single_unit_fig(2, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', example_sdf_in);
    single_unit_fig(2, cluster_i).stat_summary();
    single_unit_fig(2, cluster_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_example{cluster_i});
    single_unit_fig(2, cluster_i).set_color_options('map', color_pal.(['clu' int2str(cluster_i)]));

    % Population SDF
    single_unit_fig(3, cluster_i) = gramm('x', ops.sound_sdf_window, ...
        'y', inputSDF_plot{1}(neuron_class.(['cluster' int2str(cluster_i)]),:));
    single_unit_fig(3, cluster_i).stat_summary();
    single_unit_fig(3, cluster_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_pop{cluster_i});
    single_unit_fig(3, cluster_i).set_color_options('map', color_pal.(['clu' int2str(cluster_i)]));

    % Layout for subplots
    single_unit_fig(1, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.8 0.1 0.1], ...
        'legend', false, 'margin_height', [0 0], 'margin_width', [0 0], 'redraw', false);
    single_unit_fig(2, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.5 0.1 0.25], ...
        'margin_height', [0 0], 'margin_width', [0 0], 'redraw', false);
    single_unit_fig(3, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.2 0.1 0.25], ...
        'margin_height', [0 0], 'margin_width', [0 0], 'redraw', false);

    % Customize plot appearance
    single_unit_fig(1, cluster_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(2, cluster_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(3, cluster_i).geom_vline('xintercept', [0 413 563], 'style', {'k-', 'k--', 'k-'});
end

figure('Renderer', 'painters', 'Position', [100 100 1200 500]);
single_unit_fig.draw();

%% ===============================
%% 10. Proportion of Neurons per Area
%% ===============================

cluster_labels = {'facilitated','multiphasic','offset', 'rebound', 'suppressed', 'ramping'};
n_area = n_area(order,:);

% Normalize counts to proportions
p_n_area(:,1) = n_area(:,1)./sum(n_area(:,1));  % Auditory
p_n_area(:,2) = n_area(:,2)./sum(n_area(:,2));  % Frontal

% Plot stacked bar plot
figuren('Renderer', 'painters', 'Position', [200 200 500 350]); hold on;
bh = bar(p_n_area','stacked','LineStyle','none');
xticks([1 2]); xticklabels({'Auditory','Frontal'});
legend(cluster_labels(order), 'Location', 'eastoutside');

% Apply cluster colors
set(bh, 'FaceColor', 'Flat')
for clu_i = 1:6
    bh(clu_i).CData = color_pal.(['clu' int2str(clu_i)]);
end

%% ===============================
%% 11. Proportion analysis
%% ===============================

cluster_labels = {'facilitated','multiphasic','offset', 'rebound', 'suppressed', 'ramping'};
table(cluster_labels', n_clu_area(:,1), n_clu_area(:,2), ...
      'VariableNames', {'Cluster','Auditory','Frontal'})

for clu_i = 1:6
    contingency_table = [
        n_clu_area(clu_i,1), length(neuron_class.auditory.all) - n_clu_area(clu_i,1); 
        n_clu_area(clu_i,2), length(neuron_class.frontal.all) - n_clu_area(clu_i,2)
    ];
    [~,p] = fishertest(contingency_table);
    fprintf('Cluster %s: p = %.4f\n', cluster_labels{clu_i}, p)
end
