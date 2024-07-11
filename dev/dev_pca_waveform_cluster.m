
%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog_all, stimulusLog, spike_log] = import_exp_map();
ephysLog = clean_exp_map(ephysLog_all);

count = 0;

tic
for session_i = 1 :size(ephysLog,1)


    outfile_name = ephysLog.session{session_i}; % Processed file name
    data_in = load(fullfile(dirs.mat_data,[outfile_name '.mat']));
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), outfile_name)

    wav_names = fieldnames(data_in.spikes.waveform);

    
    for spike_i = 1:length(wav_names)
        spk_log_idx = find(strcmp(spike_log.session,outfile_name) & strcmp(spike_log.unitWAV,wav_names{spike_i}));

        if spike_log.useNeuron(spk_log_idx) == 1
            count = count + 1;
            waveform_spk_all(count,:) = nanmean(data_in.spikes.waveform.(wav_names{spike_i}));
        end
    end
end
toc

for wave_i = 1:size(waveform_spk_all,1)
    waveform_spk_norm(wave_i,:) = waveform_spk_all(wave_i,:)./max(abs(waveform_spk_all(wave_i,:)));
end


 [~, pcs, ~, ~, var_exp, ~] = pca(waveform_spk_norm);

meas=pcs;
rng('default');  % For reproducibility
eva = evalclusters(meas,'kmeans','silhouette','KList',[1:20]);

n_cl = 5; eva.OptimalK;
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

x_plot = [1 1 2 2 3];
y_plot = [1 2 1 2 1]+3;

for cl_i = 1:n_cl
    nsubplot(3,5,x_plot(cl_i),y_plot(cl_i));
    plot(nanmean(waveform_spk_norm(find(idx3==cl_i),:)),'color',[colorscale(cl_i,:) 1.0]);
    title (['n = ' int2str(length(find(idx3==cl_i)))]);
end
