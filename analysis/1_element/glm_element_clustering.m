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
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-mn',10,'-c',0.5);

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
    a = subplot(3,4,cluster_i); hold on
    plot(sdfTimes{1},nanmean(inputSDF{1}(clusterNeurons{cluster_i},:),1),'Color',plot_lineclust_color(cluster_i,:));
    vline([0 563], 'k-');  % Stimulus window
    vline([413], 'k--');   % Possible event of interest
    xlim([-150 600])
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end

%% Cluster curation
% Assign clusters to structured fields in neuron_class
neuron_class.cluster1 = [clusterNeurons{2}];
neuron_class.cluster2 = [clusterNeurons{1}];
neuron_class.cluster3 = [clusterNeurons{4}]; 
neuron_class.cluster4 = [clusterNeurons{3}];
neuron_class.cluster5 = [clusterNeurons{5}];
neuron_class.cluster6 = [clusterNeurons{6}];

% Combine all clustered neurons into a single list
neuron_class.clusterAll = [];
for cluster_i = 1:6
    neuron_class.clusterAll = [neuron_class.clusterAll ; clusterNeurons{cluster_i}];
end
neuron_class.clusterAll = sort(neuron_class.clusterAll); % Sort indices

% Identify unclustered neurons
neuron_class.unclustered = find(~ismember(1:size(sig_neurons),neuron_class.clusterAll));

% Replot average SDFs for each cluster
for cluster_i = 1:nClusters_manual
    a = subplot(3,4,cluster_i); hold on
    plot(sdfTimes{1},nanmean(inputSDF{1}(clusterNeurons{cluster_i},:),1),'Color',plot_lineclust_color(cluster_i,:));
    vline([0 563], 'k-'); 
    vline([413], 'k--'); 
    xlim([-150 600])
end

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
% Final assignment of neuron indices to neuron_class structure
neuron_class.cluster_idx.clu1 = sig_neurons(neuron_class.cluster1);
neuron_class.cluster_idx.clu2 = sig_neurons(neuron_class.cluster2);
neuron_class.cluster_idx.clu3 = sig_neurons(neuron_class.cluster3);
neuron_class.cluster_idx.clu4 = sig_neurons(neuron_class.cluster4);
neuron_class.cluster_idx.clu5 = sig_neurons(neuron_class.cluster5);
neuron_class.cluster_idx.clu6 = sig_neurons(neuron_class.cluster6);

%% 

for clu_i = 1:6

    neuron_class.cluster_idx.(['clu' int2str(clu_i) '_aud']) =...
        intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all);
    neuron_class.cluster_idx.(['clu' int2str(clu_i) '_frontal']) =...
        intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.frontal.all);

    n_clu_area(clu_i,1) = length(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.auditory.all));
    n_clu_area(clu_i,2) = length(intersect(neuron_class.cluster_idx.(['clu' int2str(clu_i)]),neuron_class.frontal.all));
end