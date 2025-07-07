% Identify significant neurons from auditory and frontal brain regions
sig_neurons = []; 
sig_neurons = [neuron_class.auditory.all; neuron_class.frontal.all];

% Clear and define input spike density function (SDF) data
clear inputSDF
inputSDF = {normFR_in.norm_fr_soundAll(sig_neurons,:)};
inputSDF_split = {normFR_in.norm_fr_soundA(sig_neurons,:),...
    normFR_in.norm_fr_soundC(sig_neurons,:),...
    normFR_in.norm_fr_soundD(sig_neurons,:),...
    normFR_in.norm_fr_soundF(sig_neurons,:),...
    normFR_in.norm_fr_soundG(sig_neurons,:)};

% Define time windows for SDF analysis
sdfTimes = {[-200:800], [-200:800], [-200:800], [-200:800], [-200:800]};
sdfEpoch = {[-50:563], [0:563], [0:563], [0:563], [0:563]};

% Define color mapping for clustering
colorMapping = [1, 1, 1, 1, 1];

% Perform consensus clustering on input SDF data
[sortIDs,idxDist, raw, respSumStruct, rawLink, myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0);
close all

% Determine the number of clusters and initialize cluster storage
nClusters_manual = myK; 
clusterNeurons = [];

% Assign neurons to their respective clusters
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end

% Plot dendrogram to visualize hierarchical clustering
figuren('Renderer', 'painters', 'Position', [100 100 500 400]);
subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
set(gca,'YTick',[],'YLim',[1 length(rawLink)]); xlabel('Similarity')

% Plot similarity matrix
subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
imagesc(raw(outPerm,outPerm));
colormap(flipud(cbrewer2('PuBu')));
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])

plot_lineclust_color = cool(nClusters_manual);

% FOR CLUSTER EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figuren('Renderer', 'painters', 'Position', [100 100 850 400]);
% Generate SDF plots for each cluster
for cluster_i = 1:nClusters_manual
    a = subplot(3,4,cluster_i); hold on
    plot(sdfTimes{1},nanmean(inputSDF{1}(clusterNeurons{cluster_i},:),1),'Color',plot_lineclust_color(cluster_i,:));
    vline([0 563], 'k-'); 
    vline([413], 'k--'); 
    xlim([-200 750])
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end

% Assign clusters to neuron_class structure
neuron_class.cluster1 = clusterNeurons{1};
neuron_class.cluster2 = clusterNeurons{2};
neuron_class.cluster3 = clusterNeurons{3};
neuron_class.cluster4 = clusterNeurons{4};
neuron_class.cluster5 = clusterNeurons{5};
neuron_class.cluster6 = clusterNeurons{6};
neuron_class.cluster7 = clusterNeurons{7};

% Generate summary SDF plots for each cluster
cluster_yaxis = {[-5 10], [-4 2], [-2 7], [-6 2], [-1 7], [-4 8], [-3 5]};
clear cluster_pop_average_sdf p_n_area n_area

for cluster_i = 1:7
    % Plot 1: Example neuron SDF plot
    cluster_pop_average_sdf(1, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', inputSDF{1}(neuron_class.(['cluster' int2str(cluster_i)]),:));
    cluster_pop_average_sdf(1, cluster_i).stat_summary();
    cluster_pop_average_sdf(1, cluster_i).axe_property('XLim', [-180 600], 'YLim', cluster_yaxis{cluster_i});
    cluster_pop_average_sdf(1, cluster_i).geom_vline('xintercept',[0 413 563]);

    n_area(cluster_i,1) = length(intersect(sig_neurons(neuron_class.(['cluster' int2str(cluster_i)])),auditory_neuron_idx));
    n_area(cluster_i,2) = length(intersect(sig_neurons(neuron_class.(['cluster' int2str(cluster_i)])),frontal_neuron_idx));

end

% Create the figure and draw the plots
figure('Renderer', 'painters', 'Position', [100 100 1250 250]);
cluster_pop_average_sdf.draw();

% Sum neuron counts per brain area
sum(n_area)

p_n_area(:,1) = n_area(:,1)./sum(n_area(:,1));
p_n_area(:,2) = n_area(:,2)./sum(n_area(:,2));

figuren('Renderer', 'painters', 'Position', [200 200 200 350]); hold on;
bar(p_n_area','stacked')
xticks([1 2]); xticklabels({'Auditory','Frontal'})
legend({'1','2','3','4','5','6','7'},'Location','southoutside')

%% 

for cluster_i = 1:4
    % Plot 1: Example neuron SDF plot

    label = ismember(sig_neurons(neuron_class.(['cluster' int2str(cluster_i)])),neuron_class.auditory.all);
    cluster_pop_average_sdf(1, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', inputSDF{1}(neuron_class.(['cluster' int2str(cluster_i)]),:), 'color', label);
    cluster_pop_average_sdf(1, cluster_i).stat_summary();
    cluster_pop_average_sdf(1, cluster_i).axe_property('XLim', [-180 600], 'YLim', cluster_yaxis{cluster_i});
    cluster_pop_average_sdf(1, cluster_i).geom_vline('xintercept',[0 413 563]);
end
figure('Renderer', 'painters', 'Position', [100 100 1000 250]);
cluster_pop_average_sdf.draw();

%% 
% Generate summary SDF plots for each cluster
cluster_yaxis = {[-2 7], [-2 7], [-6 5], [-4 2]};
clear cluster_pop_average_sdf p_n_area n_area

figuren('Renderer', 'painters', 'Position', [100 100 1000 250]);
for cluster_i = 1:4
    subplot(1,4,cluster_i); hold on
    plot(ops.sound_sdf_window, nanmean(inputSDF_split{1}(neuron_class.(['cluster' int2str(cluster_i)]),:)))
    plot(ops.sound_sdf_window, nanmean(inputSDF_split{2}(neuron_class.(['cluster' int2str(cluster_i)]),:)))
    plot(ops.sound_sdf_window, nanmean(inputSDF_split{3}(neuron_class.(['cluster' int2str(cluster_i)]),:)))
    plot(ops.sound_sdf_window, nanmean(inputSDF_split{4}(neuron_class.(['cluster' int2str(cluster_i)]),:)))
    plot(ops.sound_sdf_window, nanmean(inputSDF_split{5}(neuron_class.(['cluster' int2str(cluster_i)]),:)))
    xlim([-180 600])
    vline([0 413 563])
end

neuron_class.cluster_idx.clu1 = sort(sig_neurons(neuron_class.cluster1));
neuron_class.cluster_idx.clu2 = sort(sig_neurons(neuron_class.cluster2));
neuron_class.cluster_idx.clu3 = sort(sig_neurons(neuron_class.cluster3));
neuron_class.cluster_idx.clu4 = sort(sig_neurons(neuron_class.cluster4));
