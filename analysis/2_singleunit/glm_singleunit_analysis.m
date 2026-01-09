%% Run GLM
% This section runs a GLM-based sound modulation analysis for each neuron.
% Results are saved to disk so the analysis is not rerun unnecessarily.

% Check whether GLM results already exist
if ~exist(fullfile(dirs.root,'data','glm_modulation_data.mat'))

    % Clear variables to avoid accidental reuse
    clear glm_encoding_flag glm_beta glm_sig anova_pvalue window_sdf window_time

    % Parallel loop across neurons
    parfor neuron_i = 1:size(spike_log,1)

        % Display progress in command window
        fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1));

        % Run GLM on sound-aligned spike density function (SDF) data
        % Outputs:
        %  - glm_encoding_flag: binary flag for significant encoding
        %  - glm_beta: GLM beta weights
        %  - glm_sig: significance of GLM terms
        %  - window_sdf: SDF data used for GLM
        %  - window_time: corresponding time vector
        [glm_encoding_flag(neuron_i,:), ...
            glm_beta{neuron_i}, ...
            glm_sig{neuron_i}, ...
            window_sdf{neuron_i,1}, ...
            window_time(neuron_i,:)] = ...
            glm_sound_modulation(sdf_soundAlign_data{neuron_i});

        % Run complementary t-test based modulation analysis
        ttest_encoding_flag(neuron_i,:) = ...
            ttest_sound_modulation(sdf_soundAlign_data{neuron_i});
    end

    % Save GLM results to disk
    save(fullfile(dirs.root,'data','glm_modulation_data.mat'), ...
        'glm_encoding_flag','glm_beta','glm_sig', ...
        'window_sdf','window_time','ttest_encoding_flag','-v7.3')

else
    % Load precomputed GLM results if they already exist
    load(fullfile(dirs.root,'data','glm_modulation_data.mat'));
end


%% Identify significantly modulated neurons

% Neurons are considered significant if they encode at least one sound category
glm_sig_units = find(sum(glm_encoding_flag, 2) > 0);

% Extract the GLM time window (same for all neurons)
glm_timewin = window_time(1,:);

% Identify indices corresponding to the analysis window (0–413 ms post sound)
analysis_win_idx = find(glm_timewin >= 0 & glm_timewin <= 413);

% Total number of significantly modulated neurons
n_sig_neurons = length(glm_sig_units);


%% Classify neurons by brain region and subject

clear neuron_class

% Auditory cortex neurons with significant sound modulation
neuron_class.auditory.all = intersect(glm_sig_units, auditory_neuron_idx);

% Frontal cortex neurons with significant sound modulation
neuron_class.frontal.all  = intersect(glm_sig_units, frontal_neuron_idx);

% Neurons recorded from each monkey
neuron_class.nhp.troy = find(strcmp(spike_log.monkey,'troy'));
neuron_class.nhp.walt = find(strcmp(spike_log.monkey,'walt'));

% Count significant neurons per monkey (auditory cortex)
[ ...
    sum(strcmp(spike_log.monkey(neuron_class.auditory.all),'troy')), ...
    sum(strcmp(spike_log.monkey(neuron_class.auditory.all),'walt')) ...
    ]

% Count significant neurons per monkey (frontal cortex)
[ ...
    sum(strcmp(spike_log.monkey(neuron_class.frontal.all),'troy')), ...
    sum(strcmp(spike_log.monkey(neuron_class.frontal.all),'walt')) ...
    ]


%% Donut plots showing fraction of significant neurons

% Recompute region-specific significant neurons (redundant but explicit)
neuron_class.auditory.all = intersect(glm_sig_units, auditory_neuron_idx);
neuron_class.frontal.all  = intersect(glm_sig_units, frontal_neuron_idx);

% Create figure
figuren('Renderer','painters','Position',[680,458,169,420]);

% Auditory cortex: significant vs non-significant neurons
subplot(2,1,1)
donut([ ...
    length(neuron_class.auditory.all), ...
    length(auditory_neuron_idx) - length(neuron_class.auditory.all) ...
    ])

% Frontal cortex: significant vs non-significant neurons
subplot(2,1,2)
donut([ ...
    length(neuron_class.frontal.all), ...
    length(auditory_neuron_idx) - length(neuron_class.frontal.all) ...
    ])


%% Heatmap of normalized firing rates for significant neurons

% Load normalized firing rate data
load('normFR_in.mat')

% Compute mean firing rate in a post-sound window (200–613 ms)
metric = mean(normFR_in.norm_fr_soundAll(glm_sig_units, 200:613), 2);

% Sort neurons by response magnitude
[~, sortIdx] = sort(metric);

% Define a red–blue colormap
colorscale = abs(flipud(cbrewer2('seq','RdBu',100)));

% Plot heatmap of sound-aligned firing rates
figuren('Renderer','painters','Position',[680,458,250,420]);
imagesc(-200:800, 1:length(sortIdx), ...
    normFR_in.norm_fr_soundAll(glm_sig_units(sortIdx), :))

% Formatting
xlim([-100 800]);
ylim([1 length(sortIdx)]);
clim([-5 5])

% Mark sound onset and analysis windows
vline([0 413],'k-')     % Sound onset to analysis window end
vline([563],'k--')     % Additional reference time

% Apply colormap
colormap(colorscale)

%% GLM - identity and position
glm_identity_encode = glm_encoding_flag(:,[2:5]);
glm_position_encode = glm_encoding_flag(:,[7:10]);

glm_id_flag = sum(glm_identity_encode,2) > 0;
glm_pos_flag = sum(glm_position_encode,2) > 0;

n_aud_id_only = sum((glm_id_flag(neuron_class.auditory.all) == 1) & (glm_pos_flag(neuron_class.auditory.all) == 0));
n_aud_pos_only = sum((glm_id_flag(neuron_class.auditory.all) == 0) & (glm_pos_flag(neuron_class.auditory.all) == 1));
n_aud_both = sum((glm_id_flag(neuron_class.auditory.all) == 1) & (glm_pos_flag(neuron_class.auditory.all) == 1));
n_aud_onset = sum((glm_id_flag(neuron_class.auditory.all) == 0) & (glm_pos_flag(neuron_class.auditory.all) == 0));

n_frontal_id_only = sum((glm_id_flag(neuron_class.frontal.all) == 1) & (glm_pos_flag(neuron_class.frontal.all) == 0));
n_frontal_pos_only = sum((glm_id_flag(neuron_class.frontal.all) == 0) & (glm_pos_flag(neuron_class.frontal.all) == 1));
n_frontal_both = sum((glm_id_flag(neuron_class.frontal.all) == 1) & (glm_pos_flag(neuron_class.frontal.all) == 1));
n_frontal_onset = sum((glm_id_flag(neuron_class.frontal.all) == 0) & (glm_pos_flag(neuron_class.frontal.all) == 0));

n_aud_counts = 100*([n_aud_id_only, n_aud_pos_only, n_aud_both, n_aud_onset]./length(neuron_class.auditory.all));
n_frontal_counts = 100*([n_frontal_id_only, n_frontal_pos_only, n_frontal_both, n_frontal_onset]./length(neuron_class.frontal.all));

glm_encoding_counts = [n_aud_counts; n_frontal_counts];
figuren('Renderer', 'painters', 'Position', [100 100 400 500]); hold on;
bar([1,2], glm_encoding_counts,'stacked')
xticks([1 2])
xticklabels({'Auditory','Frontal'})
legend({'Identity','Position','Both','Onset'})


% Rows = cortex (auditory, frontal)
% Columns = neuron category (ID only, Pos only, Both, None)
glm_pos_id_cont_table = [...
    n_aud_id_only, n_aud_pos_only, n_aud_both, n_aud_onset;
    n_frontal_id_only, n_frontal_pos_only, n_frontal_both, n_frontal_onset];

% Run the chi-square test
[chi2_stat, df, p_value, expected] = chi2_independence(glm_pos_id_cont_table);

fprintf('Chi-square = %.2f, df = %d, p = %.4f\n', chi2_stat, df, p_value);
disp('Expected counts under independence:');
disp(expected);


%% GLM - number of responses
n_glm_sig_aud_id = histcounts(sum(glm_encoding_flag(neuron_class.auditory.all,[1:5]),2),-0.5:1:5.5)./length(neuron_class.auditory.all);
n_glm_sig_aud_pos = histcounts(sum(glm_encoding_flag(neuron_class.auditory.all,[6:10]),2),-0.5:1:5.5)./length(neuron_class.auditory.all);
n_glm_sig_frontal_id = histcounts(sum(glm_encoding_flag(neuron_class.frontal.all,[1:5]),2),-0.5:1:5.5)./length(neuron_class.frontal.all);
n_glm_sig_frontal_pos = histcounts(sum(glm_encoding_flag(neuron_class.frontal.all,[6:10]),2),-0.5:1:5.5)./length(neuron_class.frontal.all);

figuren('Renderer', 'painters', 'Position', [100 100 400 500]); hold on;
bar([1,2 3 4], [n_glm_sig_aud_id; n_glm_sig_aud_pos; n_glm_sig_frontal_id; n_glm_sig_frontal_pos],'stacked')
xticks([1 2 3 4 ])
xticklabels({'Id--aud','Pos--aud','Id--front','Pos--front'})
legend({'0','1','2','3','4','5'})


aud_id_counts = sum(glm_encoding_flag(neuron_class.auditory.all,[1:5]),2);
aud_pos_counts = sum(glm_encoding_flag(neuron_class.auditory.all,[6:10]),2);
front_id_counts = sum(glm_encoding_flag(neuron_class.frontal.all,[1:5]),2);
front_pos_counts = sum(glm_encoding_flag(neuron_class.frontal.all,[6:10]),2);

% Wilcoxon rank-sum test
[p_id,~,stats_id] = ranksum(aud_id_counts, front_id_counts);
[p_pos,~,stats_pos] = ranksum(aud_pos_counts, front_pos_counts);

%% Example neurons

neurons_of_interest = [789, 1214, 1261, 1303];

figuren('Renderer', 'painters', 'Position', [100 100 800 375]); hold on;

for neuron_i = 1:length(neurons_of_interest)

    clear data_in example_neuron_sdf
    neuron_idx = neurons_of_interest(neuron_i);
    data_in = sdf_soundAlign_data{neuron_idx};

    example_neuron_sdf = cell2mat(data_in(:,1));


    id_list = {'C','D','F','G'};
    pos_list = {'position_2','position_3','position_4','position_5'};

    ylim_list = {[0 40], [0 20], [0 25], [0 30]};
    
    for idx = 1:4
        clear id_trials pos_trials
        id_trials = find(strcmp(data_in(:,3),id_list{idx}));
        pos_trials = find(strcmp(data_in(:,4),pos_list{idx}));

        subplot(2,length(neurons_of_interest),neuron_i); hold on
        plot(-200:800, smooth(nanmean(example_neuron_sdf(id_trials,:)),50), 'color', [color_pal.identity(idx,:) 0.5], 'LineWidth',1.5)
        xlim([-100 600]); ylim(ylim_list{neuron_i}); axis square
        legend

        title(neuron_idx)
        subplot(2,length(neurons_of_interest),neuron_i+length(neurons_of_interest)); hold on
        plot(-200:800, smooth(nanmean(example_neuron_sdf(pos_trials,:)),50), 'color', [color_pal.position(idx,:) 0.5], 'LineWidth',1.5)
        xlim([-100 600]); ylim(ylim_list{neuron_i}); axis square
        legend
    end

end


[spike_log.session(neurons_of_interest), spike_log.unitDSP(neurons_of_interest)]