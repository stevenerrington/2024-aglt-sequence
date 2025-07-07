

parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron  
    
    sdf_in = load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data\spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
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


%% 






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

%%


%% Distance/similarity analysis
% Define the time window and extract indices for the time window from 'ops'
timewin = [-100:5:2750];  % Set the time window from -100 ms to 2750 ms
timewin_idx = find(ismember(ops.timewin, timewin));  % Get the indices of this time window from the 'ops' structure

% Define sequence elements and onset times for the experiment
seq_elements = {'A','C','D','F','G'};  % List of sequence elements (e.g., stimuli)
seq_onset_times = [0, 563, 1126, 1689, 2252];  % Onset times for each element in the sequence (in ms)

% Define sequence patterns for each sequence (set of stimuli presented in the experiment)
seq_pattern{1} = {'A','C','F','C','G'};  % Pattern for sequence 1
seq_pattern{2} = {'A','D','C','F','C'};  % Pattern for sequence 2
seq_pattern{3} = {'A','C','G','F','C'};  % Pattern for sequence 3
seq_pattern{4} = {'A','D','C','G','F'};  % Pattern for sequence 4

% Loop over 'area' cell array, which contains the names of areas (auditory, frontal)
for area = {'auditory', 'frontal'}

    % Store the spike density function (SDF) signals for each sequence in the specified time window
    signal_in = {};  % Initialize a cell array to store the signal data
    signal_in = {pca_sdf_out_seq1(neuron_class.(area{1}).all, timewin_idx),  % SDF for sequence 1 in the current area
        pca_sdf_out_seq2(neuron_class.(area{1}).all, timewin_idx),  % SDF for sequence 2
        pca_sdf_out_seq3(neuron_class.(area{1}).all, timewin_idx),  % SDF for sequence 3
        pca_sdf_out_seq4(neuron_class.(area{1}).all, timewin_idx)}; % SDF for sequence 4




    % Perform cross-condition PCA analysis to extract principal components
    [pc_out, pc_shuf_out] = get_xcond_pca(signal_in);  % Perform PCA on the input signals

    % Initialize counters and storage variables for distance calculations
    count = 0;  % Initialize count for combinations of sequences
    clear d pca_element_out_temp_a;  % Clear previous distance and temporary PCA storage

    % Loop through all pairs of sequence patterns
    for seq_i = 1:4  % Loop over first sequence (seq_i)
        for seq_j = 1:4  % Loop over second sequence (seq_j)

            count = count + 1;  % Increment the counter
            seq_comb(count,:) = [seq_i, seq_j];  % Store the current sequence pair

            % Loop over elements in the sequences
            for element_i = 1:5  % Loop over elements in the first sequence (element_i)
                for element_j = 1:5  % Loop over elements in the second sequence (element_j)

                    % Find the indices of the current sequence elements in the pattern
                    a = find(strcmp(seq_pattern{seq_i}, seq_elements{element_i}),1,'first');
                    b = find(strcmp(seq_pattern{seq_j}, seq_elements{element_j}),1,'first');

                    % Store the current element pair as a string combination
                    element_combo{element_i,element_j} = [seq_elements{element_i}, seq_elements{element_j}];

                    % If both elements are found and the sequences are different, perform Procrustes analysis
                    if ~isempty(a) && ~isempty(b) && (seq_i ~= seq_j)
                        pca_a = []; pca_b = [];  % Initialize empty arrays for PCA results
                        % Extract the PCA components for the elements of both sequences
                        element_i_idx = []; element_i_idx = find(ismember(timewin,[seq_onset_times(a):seq_onset_times(a)+413]));
                        element_j_idx = []; element_j_idx = find(ismember(timewin,[seq_onset_times(b):seq_onset_times(b)+413]));

                        % Find the size of the smallest array
                        min_length = min(length(element_i_idx), length(element_j_idx));

                        % Trim both arrays to the size of the smallest array
                        element_i_idx = element_i_idx(1:min_length);
                        element_j_idx = element_j_idx(1:min_length);


                        pca_a = pc_out{seq_i}([2,3],element_i_idx);
                        pca_b = pc_out{seq_j}([2,3],element_j_idx);

                        % Perform Procrustes analysis (measure the distance between two datasets)
                        [d(element_i,element_j,count), ~, ~] = procrustes(pca_a', pca_b');

                        % Store the PCA data for both elements in a temporary structure
                        pca_element_out_temp_a{element_i,element_j,count} = pca_a;
                        pca_element_out_temp_b{element_i,element_j,count} = pca_b;

                    else
                        % If conditions aren't met (missing elements or same sequences), set distance to NaN
                        d(element_i,element_j,count) = NaN;
                        pca_element_out_temp_a{element_i,element_j,count} = NaN;  % Set output to NaN for missing or invalid elements
                        pca_element_out_temp_b{element_i,element_j,count} = NaN;  % Ensure both entries are NaN
                    end

                end
            end
        end
    end

    % Store the results for each area (auditory or frontal) in the pca_dist structure
    pca_dist.(area{1}) = d;
    pca_dist_pc.(area{1}).a = pca_element_out_temp_a;
    pca_dist_pc.(area{1}).b = pca_element_out_temp_b;

end


% Create a figure for visualization of the PCA distance results
figuren('Renderer', 'painters', 'Position', [100 100 250 500]);

% Plot the heatmap for the auditory area (mean over the third dimension of distances)
subplot(2,1,1)
h = imagesc(tril(nanmean(pca_dist.auditory, 3)));  % Mean across the third dimension (e.g., over multiple repetitions)
clim([0 1]); title('Auditory')  % Set color limits for the heatmap
colormap(parula); box off

% Plot the heatmap for the frontal area (mean over the third dimension of distances)
subplot(2,1,2)
g = imagesc(tril(nanmean(pca_dist.frontal, 3)))  % Mean across the third dimension (e.g., over multiple repetitions)
clim([0 1]); title('Frontal')  % Set color limits for the heatmap
colormap(parula); box off


% Set transparency for zeros
alpha = ones(size(tril(nanmean(pca_dist.auditory, 3))));      % Initialize alpha (fully opaque)
alpha(tril(nanmean(pca_dist.auditory, 3)) == 0) = 0.0;        % Set transparency for zeros

% Apply the alpha data
set(h, 'AlphaData', alpha);
set(g, 'AlphaData', alpha);


%%

aud_d_tril = []; frontal_d_tril = [];
aud_d_diag = []; frontal_d_diag = [];


for seq_comb_i = [2,3,4,7,8,12]

    temp_aud = []; temp_aud = pca_dist.auditory(:,:,seq_comb_i);
    temp_frontal = []; temp_frontal = pca_dist.frontal(:,:,seq_comb_i);
    m  = tril(true(size(temp_aud)),-1);

    aud_d_tril = [aud_d_tril; temp_aud(m)];
    frontal_d_tril = [frontal_d_tril; temp_frontal(m)];

    aud_d_diag = [aud_d_diag; diag(temp_aud)];
    frontal_d_diag = [frontal_d_diag; diag(temp_frontal)];

end

pca_similiarity_data = [aud_d_tril; frontal_d_tril; aud_d_diag; frontal_d_diag];
pca_similarity_label = [repmat({'1_Aud_interelement_tril'},length(aud_d_tril),1); repmat({'3_Frontal_interelement_tril'},length(frontal_d_tril),1);...
    repmat({'2_Aud_intraelement_diag'},length(aud_d_diag),1); repmat({'4_Frontal_intraelement_diag'},length(frontal_d_diag),1)];

clear pca_similarity
pca_similarity(1,1)=gramm('x',pca_similarity_label,'y',pca_similiarity_data,'color',pca_similarity_label);
pca_similarity(1,1).stat_summary('geom',{'bar','black_errorbar'},'width',3);
fig = figure('Renderer', 'painters', 'Position', [100 100 500 400]);
pca_similarity.draw();

%%
% Store the spike density function (SDF) signals for each sequence in the specified time window
signal_in = {};  % Initialize a cell array to store the signal data
signal_in = {pca_sdf_out_seq1(neuron_class.auditory.all, timewin_idx),  % SDF for sequence 1 in the current area
    pca_sdf_out_seq2(neuron_class.auditory.all, timewin_idx),  % SDF for sequence 2
    pca_sdf_out_seq3(neuron_class.auditory.all, timewin_idx),  % SDF for sequence 3
    pca_sdf_out_seq4(neuron_class.auditory.all, timewin_idx)}; % SDF for sequence 4

% Perform cross-condition PCA analysis to extract principal components
[pc_out_auditory, pc_shuf_out_auditory] = get_xcond_pca(signal_in);  % Perform PCA on the input signals

% Store the spike density function (SDF) signals for each sequence in the specified time window
signal_in = {};  % Initialize a cell array to store the signal data
signal_in = {pca_sdf_out_seq1(neuron_class.frontal.all, timewin_idx),  % SDF for sequence 1 in the current area
    pca_sdf_out_seq2(neuron_class.frontal.all, timewin_idx),  % SDF for sequence 2
    pca_sdf_out_seq3(neuron_class.frontal.all, timewin_idx),  % SDF for sequence 3
    pca_sdf_out_seq4(neuron_class.frontal.all, timewin_idx)}; % SDF for sequence 4

% Perform cross-condition PCA analysis to extract principal components
[pc_out_frontal, pc_shuf_out_frontal] = get_xcond_pca(signal_in);  % Perform PCA on the input signals


element_1_idx = find(ismember(timewin,[0:413]));
element_2_idx = find(ismember(timewin,[563:976]));
element_3_idx = find(ismember(timewin,[1126:1539]));
element_4_idx = find(ismember(timewin,[1689:2102]));
element_5_idx = find(ismember(timewin,[2252:2665]));

figuren('Renderer', 'painters', 'Position', [100 207 634 593]);

a = subplot(2,2,1); hold on
plot(pc_out_auditory{1}(2,element_1_idx), pc_out_auditory{1}(3,element_1_idx),'linewidth',1.5,'color',[69 191 85]./255)
plot(pc_out_auditory{2}(2,element_1_idx), pc_out_auditory{2}(3,element_1_idx),'linewidth',1.5,'color',[22 128 57]./255)
plot(pc_out_auditory{3}(2,element_1_idx), pc_out_auditory{3}(3,element_1_idx),'linewidth',1.5,'color',[4 77 41]./255)
plot(pc_out_auditory{4}(2,element_1_idx), pc_out_auditory{4}(3,element_1_idx),'linewidth',1.5,'color',[0 38 28]./255)
xlim([-10 10]); ylim([-10 10]);

title('Auditory | between seq', 'FontSize', 10);

b = subplot(2,2,2); hold on
plot(pc_out_auditory{1}(2,element_1_idx), pc_out_auditory{1}(3,element_1_idx),'linewidth',1.5,'color',[114 3 79]./255)
plot(pc_out_auditory{1}(2,element_2_idx), pc_out_auditory{1}(3,element_2_idx),'linewidth',1.5,'color',[153 31 66]./255)
plot(pc_out_auditory{1}(2,element_3_idx), pc_out_auditory{1}(3,element_3_idx),'linewidth',1.5,'color',[185 63 56]./255)
plot(pc_out_auditory{1}(2,element_4_idx), pc_out_auditory{1}(3,element_4_idx),'linewidth',1.5,'color',[207 101 38]./255)
plot(pc_out_auditory{1}(2,element_5_idx), pc_out_auditory{1}(3,element_5_idx),'linewidth',1.5,'color',[215 144 0]./255)
xlim([-10 10]); ylim([-10 10]);
title('Auditory | between element', 'FontSize', 10);


c = subplot(2,2,3); hold on
plot(pc_out_frontal{1}(2,element_1_idx), pc_out_frontal{1}(3,element_1_idx),'linewidth',1.5,'color',[69 191 85]./255)
plot(pc_out_frontal{2}(2,element_1_idx), pc_out_frontal{2}(3,element_1_idx),'linewidth',1.5,'color',[22 128 57]./255)
plot(pc_out_frontal{3}(2,element_1_idx), pc_out_frontal{3}(3,element_1_idx),'linewidth',1.5,'color',[4 77 41]./255)
plot(pc_out_frontal{4}(2,element_1_idx), pc_out_frontal{4}(3,element_1_idx),'linewidth',1.5,'color',[0 38 28]./255)
xlim([-10 10]); ylim([-10 10]);
title('Frontal | between seq', 'FontSize', 10);

d = subplot(2,2,4); hold on
plot(pc_out_frontal{1}(2,element_1_idx), pc_out_frontal{1}(3,element_1_idx),'linewidth',1.5,'color',[114 3 79]./255)
plot(pc_out_frontal{1}(2,element_2_idx), pc_out_frontal{1}(3,element_2_idx),'linewidth',1.5,'color',[153 31 66]./255)
plot(pc_out_frontal{1}(2,element_3_idx), pc_out_frontal{1}(3,element_3_idx),'linewidth',1.5,'color',[185 63 56]./255)
plot(pc_out_frontal{1}(2,element_4_idx), pc_out_frontal{1}(3,element_4_idx),'linewidth',1.5,'color',[207 101 38]./255)
plot(pc_out_frontal{1}(2,element_5_idx), pc_out_frontal{1}(3,element_5_idx),'linewidth',1.5,'color',[215 144 0]./255)
xlim([-10 10]); ylim([-10 10]);
title('Frontal | between element', 'FontSize', 10);

