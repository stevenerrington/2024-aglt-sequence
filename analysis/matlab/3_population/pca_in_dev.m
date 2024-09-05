

parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron  
    
    sdf_in = load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data\spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
    event_table_in = load(fullfile(dirs.mat_data,[spike_log.session{neuron_i} '.mat']),'event_table');

    nonviol_sdf = []; nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label,'nonviol'),:);
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:,[1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:,[1000+[-200:0]])));
    
    pca_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std,50)

    for seq_i = 1:8
        pca_sdf_out_seq1(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 1,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq2(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 2,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq3(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 3,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq4(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 4,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq5(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 5,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq6(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 6,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq7(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 7,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
        pca_sdf_out_seq8(neuron_i,:) = smooth((nanmean( sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 8,:)) - baseline_fr_mean) ./ baseline_fr_std,50);
    end
end


%% Run principal components analysis


%%

pca_multiseq_sdf = [pca_sdf_out_seq1(:,1000+pca_window), pca_sdf_out_seq2(:,1000+pca_window), pca_sdf_out_seq5(:,1000+pca_window),...
    pca_sdf_out_seq6(:,1000+pca_window)];

% Marking the boundaries
% Create an index array to indicate the start of each original array
boundaryIndices = [1, cumsum(repmat(length(pca_window), 1, 4))]; % Start positions of each array in the concatenated array

clear sdf_in_regular sdf_in_shuffled 
sdf_in_regular = pca_multiseq_sdf(auditory_neuron_idx,:);

clear pcs* var*
[coeff, pcs, latent, ~, var_exp,~] = pca(sdf_in_regular');

clear pc1 pc2 pc3
pc1 = pcs(:,1); pc2 = pcs(:,2); pc3 = pcs(:,3);





clear mean_space_pos

for seq_i = 1:4
    for sound_i = 1:5
        mean_space_pos(seq_i, sound_i, 1) = nanmean(pc1(100+boundaryIndices(seq_i)+sound_times(sound_i)+[0:200]));
        mean_space_pos(seq_i, sound_i, 2) =   nanmean(pc2(100+boundaryIndices(seq_i)+sound_times(sound_i)+[0:200]));
        mean_space_pos(seq_i, sound_i, 3) =  nanmean(pc3(100+boundaryIndices(seq_i)+sound_times(sound_i)+[0:200]));

        sound_label{seq_i,sound_i} = stimulusLog.(['sound_' int2str(sound_i) '_code']){seq_i};
    end
end









% 3D PCA Plot
figuren; hold on
colorscale = flipud(cbrewer('qual','Set1',5));


for seq_i = 1:4
    for sound_i = 1:5
        switch sound_label{seq_i,sound_i}
            case 'A'
                color_in = colorscale(1,:);
            case 'C'
                color_in = colorscale(2,:);
            case 'D'
                color_in = colorscale(3,:);
            case 'F'
                color_in = colorscale(4,:);
            case 'G'
                color_in = colorscale(5,:);
        end
        scatter3(mean_space_pos(seq_i, sound_i, 1), mean_space_pos(seq_i, sound_i, 2), mean_space_pos(seq_i, sound_i, 3), 100, colorscale(sound_i,:), 'filled')

    end

            plot3(mean_space_pos(seq_i, :, 1), mean_space_pos(seq_i, :, 2), mean_space_pos(seq_i, :, 3),'k')

end

view(34.2409,7.6800)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
grid on
