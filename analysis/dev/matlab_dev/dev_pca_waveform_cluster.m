
%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- PCA waveform clusteri -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

count = 0;

tic
for session_i = 1 :size(ephysLog,1)


    outfile_name = ephysLog.session{session_i}; % Processed file name
    data_in = load(fullfile(dirs.mat_data,[outfile_name '.mat']));
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), outfile_name)

    wav_names = spike_log.unitWAV(strcmp(spike_log.session,outfile_name));

    for spike_i = 1:length(wav_names)
        spk_log_idx = find(strcmp(spike_log.session,outfile_name) & strcmp(spike_log.unitWAV,wav_names{spike_i}));

        if spike_log.useNeuron(spk_log_idx) == 1
            count = count + 1;
            waveform_in = [];
            waveform_in = nanmean(data_in.spikes.waveform.(wav_names{spike_i}))';
            
            waveform_spk_all(count,:) = waveform_in;
        end
    end
end
toc

for wave_i = 1:size(waveform_spk_all,1)
    waveform_in = [];
    waveform_in = smooth(waveform_spk_all(wave_i,:),5)';
    waveform_spk_norm(wave_i,:) = waveform_in./max(abs(waveform_in));
end


[~, pcs, ~, ~, var_exp, ~] = pca(waveform_spk_norm);

meas=pcs;
rng('default');  % For reproducibility
eva = evalclusters(meas,'kmeans','silhouette','KList',[1:20]);

n_cl = eva.OptimalK;
idx3 = kmeans(meas,n_cl,'Distance','sqeuclidean');

colorscale = cbrewer('qual','Dark2',n_cl);

figuren('Renderer', 'painters', 'Position', [100 100 1200 500]); hold on
nsubplot(3,5,[1 2 3],[1 2 3]);
for cl_i = 1:n_cl
    scatter3(pcs(find(idx3==cl_i),1),pcs(find(idx3==cl_i),2),pcs(find(idx3==cl_i),3),...
        repmat(20,sum(idx3==cl_i),1),...
        repmat(colorscale(cl_i,:),sum(idx3==cl_i),1),'filled');
end
view([-188.4000 12.9754]);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
grid on

x_plot = [1 1 2 2 3 3];
y_plot = [1 2 1 2 1 2]+3;

for cl_i = 1:n_cl
    nsubplot(3,5,x_plot(cl_i),y_plot(cl_i));
    plot(nanmean(waveform_spk_norm(find(idx3==cl_i),:)),'color',[colorscale(cl_i,:) 1.0]);
    title (['n = ' int2str(length(find(idx3==cl_i)))]);
end

ylim_clu = {[-1 1], [-1 0], [-1 1]};

figuren('Renderer', 'painters', 'Position', [100 100 1500 500]); hold on
for cl_i = 1:n_cl
    nsubplot(1,n_cl,cl_i);
    cluster_neurons = [];
    cluster_neurons = find(idx3==cl_i);

    for neu_i = 1:length(cluster_neurons)
        plot(waveform_spk_norm(cluster_neurons(neu_i),:),'color',[colorscale(cl_i,:) 0.01]);
        plot(nanmean(waveform_spk_norm(cluster_neurons,:)),'color',[colorscale(cl_i,:) 1],'LineWidth',2);
    end
end


for cl_i = 1:n_cl
    cluster_neurons = [];
    cluster_neurons = find(idx3==cl_i);

    n_clu_neuron_monkey(1,cl_i) = sum(strcmp(spike_log.monkey(cluster_neurons),'troy'));
    n_clu_neuron_monkey(2,cl_i) = sum(strcmp(spike_log.monkey(cluster_neurons),'walt'));

end