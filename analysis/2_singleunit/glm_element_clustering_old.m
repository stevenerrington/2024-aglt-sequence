% Identify significant neurons from auditory and frontal brain regions
sig_neurons = []; 
sig_neurons = [neuron_class.auditory.all; neuron_class.frontal.all]; % Combine auditory and frontal neurons

% Clear and define input spike density function (SDF) data
clear inputSDF
inputSDF = {normFR_in.norm_fr_soundAll(sig_neurons,:)}; % Get normalized firing rates for all sound trials

% Smooth SDF data for each neuron
for neuron_i = 1:size(inputSDF{1},1)
    inputSDF{1}(neuron_i,:) = smooth(inputSDF{1}(neuron_i,:),10); % Light smoothing for clustering
    inputSDF_plot{1}(neuron_i,:) = smooth(inputSDF{1}(neuron_i,:),50); % Heavier smoothing for plotting
end

% Split SDF data by different sound conditions
inputSDF_split = {normFR_in.norm_fr_soundA(sig_neurons,:),...
    normFR_in.norm_fr_soundC(sig_neurons,:),...
    normFR_in.norm_fr_soundD(sig_neurons,:),...
    normFR_in.norm_fr_soundF(sig_neurons,:),...
    normFR_in.norm_fr_soundG(sig_neurons,:)};

% Define time windows for SDF and epoch analysis
sdfTimes = {[-200:800], [-200:800], [-200:800], [-200:800], [-200:800]};
sdfEpoch = {[-50:563], [-50:563], [-50:563], [-50:563], [-50:563]};

% Define color mapping for clustering (used for each condition)
colorMapping = [1, 1, 1, 1, 1];

% Perform consensus clustering on the input SDF data
[sortIDs,idxDist, raw, respSumStruct, rawLink, myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-mn',35);

% Set number of clusters manually from clustering output
nClusters_manual = myK; 
clusterNeurons = [];

% Assign neuron indices to their respective clusters
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end

% Plot dendrogram to visualize clustering hierarchy
figuren('Renderer', 'painters', 'Position', [100 100 500 400]);
subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right'); % Draw dendrogram
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual)) % Highlight cluster boundaries
set(gca,'YTick',[],'YLim',[1 length(rawLink)]); xlabel('Similarity')

% Plot similarity matrix based on reordered dendrogram
subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic); % Mirror lower triangle to upper triangle
    end
end
imagesc(raw(outPerm,outPerm)); % Display similarity matrix
CT=cbrewer('div', 'Spectral', 100);
CT(CT > 1) = 1;

colormap(CT);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])

% Define color palette for plotting SDFs by cluster
plot_lineclust_color = cool(nClusters_manual);

% FOR CLUSTER EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figuren('Renderer', 'painters', 'Position', [100 100 850 400]);

% Plot average SDF for each cluster
for cluster_i = 1:nClusters_manual
    a = subplot(3,5,cluster_i); hold on
    plot(sdfTimes{1},nanmean(inputSDF_plot{1}(clusterNeurons{cluster_i},:),1),'Color',plot_lineclust_color(cluster_i,:));
    vline([0 563], 'k-');  % Stimulus on/off window
    vline([413], 'k--');   % Additional task-relevant time point
    xlim([-150 600])
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end

%% Cluster curation
% Assign clusters to structured fields in neuron_class (manual ordering)
neuron_class.cluster1 = [clusterNeurons{2}];
neuron_class.cluster2 = [clusterNeurons{3}];
neuron_class.cluster3 = [clusterNeurons{6}]; 
neuron_class.cluster4 = [clusterNeurons{1}];
neuron_class.cluster5 = [clusterNeurons{4}];
neuron_class.cluster6 = [clusterNeurons{5}];

% Combine all clustered neurons into a single list
neuron_class.clusterAll = [];
for cluster_i = 1:myK
    neuron_class.clusterAll = [neuron_class.clusterAll ; clusterNeurons{cluster_i}];
end
neuron_class.clusterAll = sort(neuron_class.clusterAll); % Sort indices

% Identify unclustered neurons
neuron_class.unclustered = find(~ismember(1:size(sig_neurons),neuron_class.clusterAll));

% Generate summary SDF plots using gramm for each cluster
cluster_yaxis = {[-3 3], [-5 15], [-8 5], [-4 8], [-10 4], [-5 5]};
clear cluster_pop_average_sdf p_n_area n_area

for cluster_i = 1:6
    % Use gramm to create summary plot per cluster
    cluster_pop_average_sdf(1, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', inputSDF_plot{1}(neuron_class.(['cluster' int2str(cluster_i)]),:));
    cluster_pop_average_sdf(1, cluster_i).stat_summary();
    cluster_pop_average_sdf(1, cluster_i).axe_property('XLim', [-150 600], 'YLim', cluster_yaxis{cluster_i});
    cluster_pop_average_sdf(1, cluster_i).geom_vline('xintercept',[0 413 563]); % Mark stimulus/events

    % Count neurons in auditory vs. frontal regions for each cluster
    n_area(cluster_i,1) = length(intersect(sig_neurons(neuron_class.(['cluster' int2str(cluster_i)])),auditory_neuron_idx));
    n_area(cluster_i,2) = length(intersect(sig_neurons(neuron_class.(['cluster' int2str(cluster_i)])),frontal_neuron_idx));
end

% Create figure and draw gramm plots
figure('Renderer', 'painters', 'Position', [100 100 1250 250]);
cluster_pop_average_sdf.draw();

%% 
% Final assignment of neuron indices to neuron_class structure (absolute indices)
neuron_class.cluster_idx.clu1 = sig_neurons(neuron_class.cluster1);
neuron_class.cluster_idx.clu2 = sig_neurons(neuron_class.cluster2);
neuron_class.cluster_idx.clu3 = sig_neurons(neuron_class.cluster3);
neuron_class.cluster_idx.clu4 = sig_neurons(neuron_class.cluster4);
neuron_class.cluster_idx.clu5 = sig_neurons(neuron_class.cluster5);
neuron_class.cluster_idx.clu6 = sig_neurons(neuron_class.cluster6);

%% 
clustered_neurons = [];

% Compute auditory / frontal neuron membership per cluster
for clu_i = 1:6

    neuron_class.cluster_idx.(['clu' int2str(clu_i) '_aud']) =...
        intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),auditory_neuron_idx);
    neuron_class.cluster_idx.(['clu' int2str(clu_i) '_frontal']) =...
        intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),frontal_neuron_idx);

    n_clu_area(clu_i,1) = length(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),auditory_neuron_idx));
    n_clu_area(clu_i,2) = length(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),frontal_neuron_idx));

    clustered_neurons = [clustered_neurons; neuron_class.cluster_idx.(['clu' int2str(clu_i)])];
end

% Normalize counts relative to total auditory/frontal pools
p_clu_area(:,1) = n_clu_area(:,1)./length(neuron_class.auditory.all);
p_clu_area(:,2) = n_clu_area(:,2)./length(neuron_class.frontal.all);

%% 
% Identify neurons that did not fall into any cluster
clustered_neurons = sort(clustered_neurons);
unclustered_neurons = sig_neurons(find(~ismember(sig_neurons, clustered_neurons)));

% Count unclustered auditory / frontal neurons
length(intersect(unclustered_neurons,auditory_neuron_idx))
length(intersect(unclustered_neurons,frontal_neuron_idx))

% Define indices of example neurons for each cluster
clu_exp_neurons = [18 9 8 12 2 30];

% Set y-axis limits for the example neuron plots (SDF)
plot_ylim_example = {[2 8], [0 10], [0 50], [0 15], [6 14], [30 50]};

% Set y-axis limits for the population SDF plots
plot_ylim_pop = {[-4 4], [-5 20], [-10 5], [-2 8], [-10 4], [-3 6]};

% Define horizontal positioning for each cluster plot in the figure layout
plot_startxpos = [0.1 0.25 0.4 0.55 0.7 0.85];

% Define the desired order of clusters for plotting
order = [1 2 3 4 5 6];

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
cluster_labels = {'facilitated','multiphasic','offset', 'rebound', 'suppressed', 'ramping'};

n_area = n_area(order,:);
% Calculate total neuron count per brain area
sum(n_area)

% Normalize neuron counts to get proportions for each area
p_n_area(:,1) = n_area(:,1)./sum(n_area(:,1));  % Auditory
p_n_area(:,2) = n_area(:,2)./sum(n_area(:,2));  % Frontal

% Create bar plot showing proportion of each cluster per area
figuren('Renderer', 'painters', 'Position', [200 200 500 350]); hold on;
bh = bar(p_n_area','stacked','LineStyle','none');  % Transpose to plot area as x-axis
xticks([1 2]); xticklabels({'Auditory','Frontal'})  % Label x-axis
legend(cluster_labels(order),'Location','eastoutside')  % Add legend for cluster types

% Apply cluster-specific colors to the bars
set(bh, 'FaceColor', 'Flat')
for clu_i = 1:6
    bh(clu_i).CData = color_pal.(['clu' int2str(clu_i)]);
end
