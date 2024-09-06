
parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1))

    sound_list = {'A','C','D','F','G'};
    for sound_i = 1:5
        sdf_in = []; sdf_in = cell2mat(sdf_soundAlign_data{neuron_i}(:,1))
        baseline_in = []; baseline_in = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),sound_list{sound_i}),201+[-200:0]),2);
        activity_in = []; activity_in = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),sound_list{sound_i}),200+[0:400]),2);

        test(neuron_i,sound_i) = ttest(baseline_in, activity_in);

    end

end

sig_neurons = []; sig_neurons = find(sum(test,2) > 0);

inputSDF = {norm_fr_soundA(sig_neurons,:), norm_fr_soundC(sig_neurons,:), norm_fr_soundD(sig_neurons,:), norm_fr_soundF(sig_neurons,:), norm_fr_soundG(sig_neurons,:)};

sdfTimes = {[-200:600], [-200:600], [-200:600], [-200:600], [-200:600]};
sdfEpoch = {[0:400], [0:400], [0:400], [0:400], [0:400]};

colorMapping = [1, 1, 1, 1, 1];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0.5);

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
colormap(flipud(bone));
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])


% Generate a quick sdf of each cluster
for cluster_i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    a = subplot(1,1,1); hold on
    plot(sdfTimes{1},nanmean(inputSDF{1}(clusterNeurons{cluster_i},:),1));
    plot(sdfTimes{1},nanmean(inputSDF{2}(clusterNeurons{cluster_i},:),1));
    plot(sdfTimes{1},nanmean(inputSDF{3}(clusterNeurons{cluster_i},:),1));
    plot(sdfTimes{1},nanmean(inputSDF{4}(clusterNeurons{cluster_i},:),1));
    plot(sdfTimes{1},nanmean(inputSDF{5}(clusterNeurons{cluster_i},:),1));
    vline([0], 'k--'); xlim([-200 600])


    xlim([-200 600])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end
















sig_list(:,1) = (sum(test,2) > 0);
sig_list(:,2) = sum(sum(encoding_flag,2) > 0);




joint_sig = intersect(find((sum(test,2) > 0)),find((sum(encoding_flag,2) > 0)))

clear glm_output encoding_flag encoding_beta sdf_normalized

parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1))

    [glm_output{neuron_i}, encoding_flag(neuron_i,:), encoding_beta(neuron_i,:), sdf_normalized{neuron_i,1}] =...
        glm_sound_modulation(sdf_soundAlign_data{neuron_i});
end

for neuron_i = 1:size(spike_log,1)
    mean_fr = sdf_normalized{neuron_i,1}

    norm_fr_baseline(neuron_i,:) = nanmean(sdf_normalized{neuron_i,1}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),:));
    norm_fr_soundA(neuron_i,:) = nanmean(sdf_normalized{neuron_i,1}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'),:));
    norm_fr_soundC(neuron_i,:) = nanmean(sdf_normalized{neuron_i,1}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'),:));
    norm_fr_soundD(neuron_i,:) = nanmean(sdf_normalized{neuron_i,1}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'),:));
    norm_fr_soundF(neuron_i,:) = nanmean(sdf_normalized{neuron_i,1}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'),:));
    norm_fr_soundG(neuron_i,:) = nanmean(sdf_normalized{neuron_i,1}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'),:));
end



neurons_in = find((sum(test,2) > 0));
auditory_neuron_in = intersect(neurons_in, auditory_neuron_idx);


for neuron_i = 1:10
figuren; hold on
plot(ops.sound_sdf_window,smooth(norm_fr_baseline(joint_sig(neuron_i),:),50),'k','LineWidth',2)
plot(ops.sound_sdf_window,smooth(norm_fr_soundA(joint_sig(neuron_i),:),50))
plot(ops.sound_sdf_window,smooth(norm_fr_soundC(joint_sig(neuron_i),:),50))
plot(ops.sound_sdf_window,smooth(norm_fr_soundD(joint_sig(neuron_i),:),50))
plot(ops.sound_sdf_window,smooth(norm_fr_soundF(joint_sig(neuron_i),:),50))
plot(ops.sound_sdf_window,smooth(norm_fr_soundG(joint_sig(neuron_i),:),50))
end




figuren;
subplot(1,2,1); hold on
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundA(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundC(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundD(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundF(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundG(intersect(neurons_in,frontal_neuron_idx),:)),50))

subplot(1,2,2); hold on
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundA(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundC(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundD(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundF(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundG(intersect(neurons_in,auditory_neuron_idx),:)),50))



length(intersect(neurons_in,frontal_neuron_idx))
length(intersect(neurons_in,auditory_neuron_idx))

