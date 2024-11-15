% Loop through each neuron
parfor neuron_i = 1:size(spike_log,1)
    % Display progress for the current neuron
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); 
    
    % Load the spike data for the current neuron
    sdf_in = load(fullfile(dirs.root,'data','spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
    
    % Load the event table for the current session
    event_table_in = load(fullfile(dirs.mat_data, [spike_log.session{neuron_i} '.mat']), 'event_table');

    % Extract spike density function (SDF) for 'nonviolent' condition
    nonviol_sdf = []; 
    nonviol_sdf = sdf_in.sdf.sequenceOnset(strcmp(event_table_in.event_table.cond_label, 'nonviol'), :);
    
    % Calculate baseline firing rate mean and standard deviation over a pre-stimulus period (-200 ms to 0 ms)
    baseline_fr_mean = nanmean(nanmean(nonviol_sdf(:, [1000+[-200:0]])));
    baseline_fr_std = nanstd(nanmean(nonviol_sdf(:, [1000+[-200:0]])));

    % Compute and smooth normalized SDF for the current neuron across time points
    pca_sdf_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mean) ./ baseline_fr_std, 50);

    % Calculate and smooth the SDF for each sequence condition, normalized by baseline firing rate
    pca_sdf_out_seq1(neuron_i,:) = smooth((nanmean(sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 1 | event_table_in.event_table.cond_value == 5, :)) - baseline_fr_mean) ./ baseline_fr_std, 50);
    pca_sdf_out_seq2(neuron_i,:) = smooth((nanmean(sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 2 | event_table_in.event_table.cond_value == 6, :)) - baseline_fr_mean) ./ baseline_fr_std, 50);
    pca_sdf_out_seq3(neuron_i,:) = smooth((nanmean(sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 3 | event_table_in.event_table.cond_value == 7, :)) - baseline_fr_mean) ./ baseline_fr_std, 50);
    pca_sdf_out_seq4(neuron_i,:) = smooth((nanmean(sdf_in.sdf.sequenceOnset(event_table_in.event_table.cond_value == 4 | event_table_in.event_table.cond_value == 8, :)) - baseline_fr_mean) ./ baseline_fr_std, 50);

end



%% All sequences
pc_out_auditory = perform_pca_and_plot(neuron_class.auditory.all, pca_sdf_out);
pc_out_frontal = perform_pca_and_plot(neuron_class.frontal.all, pca_sdf_out);

peri_element_timewin = [0:413]; 
figuren('Renderer', 'painters', 'Position', [100 100 900 400]); 
subplot(1,2,1); hold on
plot(pc_out_auditory.obs.pcs(100+0+peri_element_timewin,2),pc_out_auditory.obs.pcs(100+0+peri_element_timewin,3))
plot(pc_out_auditory.obs.pcs(100+563+peri_element_timewin,2),pc_out_auditory.obs.pcs(100+563+peri_element_timewin,3))
plot(pc_out_auditory.obs.pcs(100+1126+peri_element_timewin,2),pc_out_auditory.obs.pcs(100+1126+peri_element_timewin,3))
plot(pc_out_auditory.obs.pcs(100+1689+peri_element_timewin,2),pc_out_auditory.obs.pcs(100+1689+peri_element_timewin,3))
plot(pc_out_auditory.obs.pcs(100+2252+peri_element_timewin,2),pc_out_auditory.obs.pcs(100+2252+peri_element_timewin,3))
xlabel('PC2'); ylabel('PC3')

subplot(1,2,2); hold on
plot(pc_out_frontal.obs.pcs(100+0+peri_element_timewin,1),pc_out_frontal.obs.pcs(100+0+peri_element_timewin,2))
plot(pc_out_frontal.obs.pcs(100+563+peri_element_timewin,1),pc_out_frontal.obs.pcs(100+563+peri_element_timewin,2))
plot(pc_out_frontal.obs.pcs(100+1126+peri_element_timewin,1),pc_out_frontal.obs.pcs(100+1126+peri_element_timewin,2))
plot(pc_out_frontal.obs.pcs(100+1689+peri_element_timewin,1),pc_out_frontal.obs.pcs(100+1689+peri_element_timewin,2))
plot(pc_out_frontal.obs.pcs(100+2252+peri_element_timewin,1),pc_out_frontal.obs.pcs(100+2252+peri_element_timewin,2))
xlabel('PC1'); ylabel('PC2')


%% 
timewin = [-100:1:2750];
timewin_idx = find(ismember(ops.timewin, timewin));  % Get the index of the time window

% Store SDF signals for each sequence in the specified time window
signal_in = {};
signal_in = {pca_sdf_out_seq1(neuron_class.auditory.all, timewin_idx), 
             pca_sdf_out_seq2(neuron_class.auditory.all, timewin_idx), 
             pca_sdf_out_seq3(neuron_class.auditory.all, timewin_idx), 
             pca_sdf_out_seq4(neuron_class.auditory.all, timewin_idx)};


% Perform cross-condition PCA analysis
[pc_out, pc_shuf_out] = get_xcond_pca(signal_in);


pc_idx_plot = 1;

figuren;
plot(timewin, pc_out{1}(pc_idx_plot,:))
plot(timewin, pc_out{2}(pc_idx_plot,:))
plot(timewin, pc_out{3}(pc_idx_plot,:))
plot(timewin, pc_out{4}(pc_idx_plot,:))
vline(sound_onset_ms,'k')
vline(sound_onset_ms+413,'k--')

for seq_i = 1:4
    for seq_j = 1:4
        pca_corr_matrix(seq_i, seq_j) = corr(pc_out{seq_i}(1,:)',pc_out{seq_j}(1,:)')
    end
end

figure
heatmap(pca_corr_matrix)


%% Generate subsamples from the population
n_obs_psuedopop = 1000;  % Number of observations in the pseudopopulation
subsamples = gen_psudeo_population(neuron_class.auditory.all, n_obs_psuedopop);  % Generate pseudopopulation subsamples

% Define time window for analysis (-200 ms to 3000 ms)
timewin = [-100:1:2750];
timewin_idx = find(ismember(ops.timewin, timewin));  % Get the index of the time window

% Initialize the first subsample
subsample_i = 1;

% Store SDF signals for each sequence in the specified time window
signal_in = {};
signal_in = {pca_sdf_out_seq1(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq2(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq3(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq4(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq5(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq6(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq7(subsamples{subsample_i}, timewin_idx), 
             pca_sdf_out_seq8(subsamples{subsample_i}, timewin_idx)};


% Perform cross-condition PCA analysis
[pc_out, pc_shuf_out] = get_xcond_pca(signal_in);


n_ele = 5;
seq_i = 1;

for seq_i = 1:4
    pca_co_seq{seq_i}(1,:) = [nanmean(pc_out{seq_i}(1,100+0+peri_element_timewin)), nanmean(pc_out{seq_i}(2,100+0+peri_element_timewin)), nanmean(pc_out{seq_i}(3,100+0+peri_element_timewin))];
    pca_co_seq{seq_i}(2,:) = [nanmean(pc_out{seq_i}(1,100+563+peri_element_timewin)), nanmean(pc_out{seq_i}(2,100+563+peri_element_timewin)), nanmean(pc_out{seq_i}(3,100+563+peri_element_timewin))];
    pca_co_seq{seq_i}(3,:) = [nanmean(pc_out{seq_i}(1,100+1126+peri_element_timewin)), nanmean(pc_out{seq_i}(2,100+1126+peri_element_timewin)), nanmean(pc_out{seq_i}(3,100+1126+peri_element_timewin))];
    pca_co_seq{seq_i}(4,:) = [nanmean(pc_out{seq_i}(1,100+1689+peri_element_timewin)), nanmean(pc_out{seq_i}(2,100+1689+peri_element_timewin)), nanmean(pc_out{seq_i}(3,100+1689+peri_element_timewin))];
    pca_co_seq{seq_i}(5,:) = [nanmean(pc_out{seq_i}(1,100+2252+peri_element_timewin)), nanmean(pc_out{seq_i}(2,100+2252+peri_element_timewin)), nanmean(pc_out{seq_i}(3,100+2252+peri_element_timewin))];
end

colormap_in = parula(5);

figuren; hold on

sound_labels = {'A','C','D','F','G'};

for seq_i = 1:4
    scatter3(pca_co_seq{seq_i}(1,1),pca_co_seq{seq_i}(1,2),pca_co_seq{seq_i}(1,3),100,'filled','MarkerFaceColor', colormap_in(find(strcmp(stimulusLog.sound_1_code{seq_i},sound_labels)),:))
    scatter3(pca_co_seq{seq_i}(2,1),pca_co_seq{seq_i}(2,2),pca_co_seq{seq_i}(2,3),100,'filled','MarkerFaceColor', colormap_in(find(strcmp(stimulusLog.sound_2_code{seq_i},sound_labels)),:))
    scatter3(pca_co_seq{seq_i}(3,1),pca_co_seq{seq_i}(3,2),pca_co_seq{seq_i}(3,3),100,'filled','MarkerFaceColor', colormap_in(find(strcmp(stimulusLog.sound_3_code{seq_i},sound_labels)),:))
    scatter3(pca_co_seq{seq_i}(4,1),pca_co_seq{seq_i}(4,2),pca_co_seq{seq_i}(4,3),100,'filled','MarkerFaceColor', colormap_in(find(strcmp(stimulusLog.sound_4_code{seq_i},sound_labels)),:))
    scatter3(pca_co_seq{seq_i}(5,1),pca_co_seq{seq_i}(5,2),pca_co_seq{seq_i}(5,3),100,'filled','MarkerFaceColor', colormap_in(find(strcmp(stimulusLog.sound_5_code{seq_i},sound_labels)),:))

end

xlabel('PC1'); xlim([-50 50])
ylabel('PC2'); ylim([-50 50])
zlabel('PC3'); zlim([-50 50])