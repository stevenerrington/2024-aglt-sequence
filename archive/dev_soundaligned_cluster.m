
clear auditory_neuron_idx inputSDF
auditory_neuron_idx = find(strcmp(spike_log.area,'R') | strcmp(spike_log.area,'A1'));
inputSDF = {norm_fr_soundA(auditory_neuron_idx,:)};

sdfTimes = {[-200:600]};
sdfEpoch = {[0:400]};

colorMapping = [1];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0.5);

%%
close all
nClusters_manual = myK; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end


% Plot dendrogram
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
set(gca,'YTick',[],'YLim',[1 length(rawLink)]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
imagesc(raw(outPerm,outPerm));
colormap(flipud(summer));
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])


% Generate a quick sdf of each cluster
for cluster_i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    a = subplot(1,1,1); hold on
    plot(sdfTimes{1},nanmean(norm_fr_soundA(clusterNeurons{cluster_i},:),1), 'color', 'k');
    vline([0], 'k--'); xlim([-200 600])


    xlim([-200 600])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end

cluster_i = 8;

neurons_in = clusterNeurons{cluster_i};

for neuron_i = 1:length(neurons_in)
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    a = subplot(1,1,1); hold on
    plot(sdfTimes{1},norm_fr_soundA(neurons_in(neuron_i),:), 'color', 'k');
    vline([0], 'k--'); xlim([-200 600])

    xlim([-200 600])
    
    title(['Neuron ' int2str(neuron_i) ])

end

figuren;
plot(sdfTimes{1},nanmean(norm_fr_soundA(neurons_in,:)), 'color', 'k');
