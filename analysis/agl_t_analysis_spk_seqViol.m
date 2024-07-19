

%% IN DEV: 2024-07-08
% SE-  just developed some code to try and form an alignment point for the
% time of violation. Currently just single violation, and this is collapsed
% across all conditions


%% Setup workspace
% Set parameters
ops.timewin = -1000:5000;
ops.bl_win = -150:-50;
ops.freq = [2 200];
ops.sdf_filter = 'PSP';
ops.sig_threshold = 2;
ops.min_sig_time = 50;
ops.freq = [2 60];

% Set parameters for sound specific stimuli
ops.sound_sdf_window = -200:600;
sound_onset_ms = [0, 563, 1126, 1689, 2252];
zero_offset = abs(ops.timewin(1));


%% Extract sound aligned SDF
% Load data
neuron_count = 0;
lfp_count = 0;

for session_i = 1:size(ephysLog,1)
    datafile = ephysLog.session{session_i};
    load(fullfile(dirs.mat_data,datafile))
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), datafile)

    % Align spikes and generate SDF
    aligntime = event_table.stimulusOnset_ms;
    [sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);
    [~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);

    % Find all neurons in given session
    neuron_list = spike_log.unitDSP(strcmp(spike_log.session,datafile));

    for ch_i = 1:32
        sound_lfp = {};
        count= 0;
        lfp_count = lfp_count + 1;

        for trial_i = 1:size(event_table,1)
            if ~strcmp(event_table.cond_label(trial_i),'error')
                for sound_i = 1:5
                    count = count + 1;
                    sound_lfp{count,1} = lfp_aligned(ch_i,zero_offset+ops.sound_sdf_window+sound_onset_ms(sound_i),trial_i);
                    sound_lfp{count,2} = stimulusLog.(['sound_' int2str(sound_i) '_code']){event_table.cond_value(trial_i)};
                    sound_lfp{count,3} = ['position_' int2str(sound_i)];
                    sound_lfp{count,4} = event_table.cond_label(trial_i);
                end
            end
        end

        sound_lfp_out{lfp_count} = sound_lfp;
    end

    % For each neuron, get the sound aligned SDF, and save relevant information
    % in the sound_sdf variable
    for neuron_i = 1:length(neuron_list)
        neuron_label = neuron_list{neuron_i};
        neuron_count = neuron_count + 1;

        sound_sdf = {};
        count= 0;

        for trial_i = 1:size(event_table,1)
            if ~strcmp(event_table.cond_label(trial_i),'error')
                for sound_i = 1:5
                    count = count + 1;
                    sound_sdf{count,1} = sdf.(neuron_label)(trial_i,zero_offset+ops.sound_sdf_window+sound_onset_ms(sound_i));
                    sound_sdf{count,2} = raster.(neuron_label){trial_i}-sound_onset_ms(sound_i);
                    sound_sdf{count,3} = stimulusLog.(['sound_' int2str(sound_i) '_code']){event_table.cond_value(trial_i)};
                    sound_sdf{count,4} = ['position_' int2str(sound_i)];
                    sound_sdf{count,5} = event_table.cond_label(trial_i);
                end
            end
        end

        sound_sdf_out{neuron_count} = sound_sdf;
        neuron_baseline(neuron_count,:) = nanmean(sdf.(neuron_label)(:,zero_offset+ops.bl_win));

    end
end


%% Normalize A-sound aligned spike density function
for neuron_i = 1:size(spike_log,1)
    a_trials = strcmp(sound_sdf_out{neuron_i}(:,3),'A');
    c_trials = strcmp(sound_sdf_out{neuron_i}(:,3),'C');
    g_trials = strcmp(sound_sdf_out{neuron_i}(:,3),'G');
    f_trials = strcmp(sound_sdf_out{neuron_i}(:,3),'F');
    d_trials = strcmp(sound_sdf_out{neuron_i}(:,3),'D');

    baseline_fr_mu = nanmean(neuron_baseline(neuron_i,:));
    baseline_fr_std = nanstd(neuron_baseline(neuron_i,:));

    norm_fr_soundA(neuron_i,:) = (nanmean(cell2mat(sound_sdf_out{neuron_i}(a_trials,1)))-baseline_fr_mu)./baseline_fr_std;
    norm_fr_soundC(neuron_i,:) = (nanmean(cell2mat(sound_sdf_out{neuron_i}(c_trials,1)))-baseline_fr_mu)./baseline_fr_std;
    norm_fr_soundG(neuron_i,:) = (nanmean(cell2mat(sound_sdf_out{neuron_i}(g_trials,1)))-baseline_fr_mu)./baseline_fr_std;
    norm_fr_soundF(neuron_i,:) = (nanmean(cell2mat(sound_sdf_out{neuron_i}(f_trials,1)))-baseline_fr_mu)./baseline_fr_std;
    norm_fr_soundD(neuron_i,:) = (nanmean(cell2mat(sound_sdf_out{neuron_i}(d_trials,1)))-baseline_fr_mu)./baseline_fr_std;

    norm_fr_sound_all(neuron_i,:) = (nanmean(cell2mat(sound_sdf_out{neuron_i}(:,1)))-baseline_fr_mu)./baseline_fr_std;
end



%%






























% 
% 
% 
% 
% 
% 
% 
% %% Run principal components analysis
% auditory_neuron_idx = find(strcmp(spike_log.area,'R') | strcmp(spike_log.area,'A1'));
% frontal_neuron_idx = find(strcmp(spike_log.area,'44') | strcmp(spike_log.area,'45'));
% 
% sdf_in_regular = norm_fr_sound_all(auditory_neuron_idx,:);
% 
% for neuron_i = 1:size(sdf_in_regular,1)
%     sdf_in_shuffled(neuron_i,:) = sdf_in_regular(neuron_i,randperm(size(sdf_in_regular,2)));
% end
% 
% clear pcs* var*
% [coeff, pcs, latent, ~, var_exp,~] = pca(sdf_in_regular');
% [~, pcs_shuffled, ~, ~, var_exp_shuffled,~] = pca(sdf_in_shuffled');
% 
% pca_color = [162 13 30; 172 40 55; 183 67 80; 193 94 105; 203 121 130; 214 147 155; 224 174 180; 234 201 205]./255;
% pca_color_shuff = [140 140 140; 153 153 153; 166 166 166; 178 178 178; 191 191 191; 204 204 204; 217 217 217]./255;
% 
% 
% [var_exp(1:8)';var_exp_shuffled(1:8)']
% 
% clear pc1 pc2 pc3
% pc1 = pcs(:,1); pc2 = pcs(:,2); pc3 = pcs(:,3);
% pc4 = pcs(:,4); pc5 = pcs(:,5); pc6 = pcs(:,6);
% onset_time_idx = abs(ops.sound_sdf_window(1));
% offset_time_idx = abs(ops.sound_sdf_window(1))+413;
% onset2_time_idx = abs(ops.sound_sdf_window(1))+563;
% 
% 
% figuren('Renderer', 'painters', 'Position', [100 100 1300 400]); hold on; 
% % Variance explained
% n_vars = 6;
% colorscale = flipud(cbrewer('seq','PuRd',n_vars));
% colorscale_shuf = flipud(cbrewer('seq','Greys',n_vars));
% 
% 
% nsubplot(3,10,[1 2 3],[1 2]); hold on
% b_obs = bar([1:n_vars],var_exp([1:n_vars]),'LineStyle','None','FaceAlpha',0.5);
% b_shuf = bar([1:n_vars],var_exp_shuffled([1:n_vars]),'LineStyle','None','FaceAlpha',0.5);
% xlabel('Principal Component')
% ylabel('Cumulative Variance Explained (%)')
% xticks([1:1:n_vars])
% ylim([0 50])
% 
% for var_i = 1:n_vars
%     b_obs.FaceColor = 'flat';
%     b_shuf.FaceColor = 'flat';
%     b_obs.CData(var_i,:) = colorscale(var_i,:);
%     b_shuf.CData(var_i,:) = colorscale_shuf(1,:);
% end
% 
% 
% % PCA x time
% nsubplot(3,10,[1], [4 5 6]); hold on
% plot(ops.sound_sdf_window, pc1, 'color',colorscale(1,:),'LineWidth',1.5)
% vline([0 413 563],'k')
% set(gca,'Xcolor',[1 1 1])
% 
% nsubplot(3,10,[2], [4 5 6]); hold on
% plot(ops.sound_sdf_window, pc2, 'color',colorscale(2,:),'LineWidth',1.5)
% vline([0 413 563],'k'); ylabel('PC')
% set(gca,'Xcolor',[1 1 1])
% 
% nsubplot(3,10,[3], [4 5 6]); hold on
% plot(ops.sound_sdf_window, pc3, 'color',colorscale(3,:),'LineWidth',1.5)
% xlabel('Time from stimulus onset (ms)')
% vline([0 413 563],'k')
% 
% % 3D PCA Plot
% nsubplot(3,10,[1 2 3], [8 9 10]); hold on 
% color_line3(pc1, pc2, pc3, ops.sound_sdf_window,'LineWidth',2)
% scatter3(pc1(onset_time_idx),pc2(onset_time_idx),pc3(onset_time_idx),75,'k','^','filled')
% scatter3(pc1(offset_time_idx),pc2(offset_time_idx),pc3(offset_time_idx),75,'k','v','filled')
% scatter3(pc1(onset2_time_idx),pc2(onset2_time_idx),pc3(onset2_time_idx),75,'k','o','filled')
% view(34.2409,7.6800)
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); 
% xlim([-40 40]); ylim([-30 30]); zlim([-30 30])
% 
% 
% grid on
% 
% % /////////////////////////////////////
% figuren('Renderer', 'painters', 'Position', [100 100 500 900]); hold on; 
% 
% % PCA x time
% nsubplot(4,1,[1], [1]); hold on
% plot(ops.sound_sdf_window, pc4, 'color',colorscale(1,:),'LineWidth',1.5)
% vline([0 413 563],'k')
% set(gca,'Xcolor',[1 1 1])
% 
% nsubplot(4,1,[2], [1]); hold on
% plot(ops.sound_sdf_window, pc5, 'color',colorscale(2,:),'LineWidth',1.5)
% vline([0 413 563],'k'); ylabel('PC')
% set(gca,'Xcolor',[1 1 1])
% 
% nsubplot(4,1,[3], [1]); hold on
% plot(ops.sound_sdf_window, pc6, 'color',colorscale(3,:),'LineWidth',1.5)
% vline([0 413 563],'k'); ylabel('PC')
% set(gca,'Xcolor',[1 1 1])
% 
% nsubplot(4,1,[4], [1]); hold on
% plot(ops.sound_sdf_window, pcs(:,7), 'color',colorscale(4,:),'LineWidth',1.5)
% xlabel('Time from stimulus onset (ms)')
% vline([0 413 563],'k')
% 
% % /////////////////////////////////////
% figuren('Renderer', 'painters', 'Position', [100 100 500 400]); hold on; 
% 
% % PCA x time
% nsubplot(1,1,[1], [1]); hold on
% plot(ops.sound_sdf_window,  pcs(:,8), 'color',colorscale(1,:),'LineWidth',1.5)
% vline([0 413 563],'k')
% 
% 
% %%
% [coeff, pcs, latent, ~, var_exp,~] = pca(sdf_in_regular);
% 
% meas=pcs;
% rng('default');  % For reproducibility
% eva = evalclusters(meas,'kmeans','silhouette','KList',[1:20]);
% 
% n_cl = 3; eva.OptimalK;
% idx3 = kmeans(meas,n_cl,'Distance','sqeuclidean');
% 
% colorscale = jet(n_cl);
% 
% figuren('Renderer', 'painters', 'Position', [100 100 1200 500]); hold on
% nsubplot(3,5,[1 2 3],[1 2 3]);
% for cl_i = 1:n_cl
%     scatter3(pcs(find(idx3==cl_i),1),pcs(find(idx3==cl_i),2),pcs(find(idx3==cl_i),3),...
%         repmat(20,sum(idx3==cl_i),1),...
%         repmat(colorscale(cl_i,:),sum(idx3==cl_i),1),'filled');
% end
% view([-188.4000 12.9754]);
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
% grid on
% 
% x_plot = [1 1 2 2 3];
% y_plot = [1 2 1 2 1]+3;
% 
% figuren('Renderer', 'painters', 'Position', [100 100 1400 400]); hold on;
% for cl_i = 1:n_cl
%     nsubplot(1,n_cl,1,cl_i);
%     plot(ops.sound_sdf_window, nanmean(sdf_in_regular(find(idx3==cl_i),:)),'color','k');
%     title (['n = ' int2str(length(find(idx3==cl_i)))]);
%     vline([0 413 563],'r')
% end
% 
% clu_i = 1;
% cluster_neurons = find(idx3==clu_i);
% 
% for cl_neur = 1:5 %length(cluster_neurons)
%     figuren;
%     plot(ops.sound_sdf_window, smooth(sdf_in_regular(cluster_neurons(cl_neur),:),10),'color','k');
% end
% 
% 
% 
% 
% %%
% % Determine periods of significant modulation
% sig_onset_times = get_sig_onset_sound(ops, sound_sdf_out{neuron_i}, neuron_baseline(neuron_i,:));
