
%% Run GLM
% Parallel loop to perform GLM analysis for each neuron
if ~exist(fullfile(dirs.root,'data','glm_modulation_data.mat'))
    clear glm_encoding_flag
    parfor neuron_i = 1:size(spike_log,1)
        fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

        % Perform GLM analysis on sound-aligned SDF data for the current neuron
        [glm_encoding_flag(neuron_i,:), glm_beta{neuron_i}, glm_sig{neuron_i}, anova_pvalue{neuron_i},...
            window_sdf{neuron_i,1}, window_time(neuron_i,:)] = glm_sound_modulation(sdf_soundAlign_data{neuron_i});
        [ttest_encoding_flag(neuron_i,:)] = ttest_sound_modulation(sdf_soundAlign_data{neuron_i});
    end
    save(fullfile(dirs.root,'data','glm_modulation_data.mat'),'glm_encoding_flag','glm_beta','glm_sig','anova_pvalue','window_sdf','ttest_encoding_flag','-v7.3')
else
    load(fullfile(dirs.root,'data','glm_modulation_data.mat'));
end

% Identify neurons with significant modulation based on GLM results
glm_sig_units = find(sum(glm_encoding_flag, 2) > 0); % Find neurons with significant encoding in any sound category

% For future use: state the time windows for the GLM
glm_timewin = -150:10:750; % Define the time window for plotting beta weights
analysis_win_idx = find(glm_timewin >= 0 & glm_timewin <= 413); % Find the relevant indicies for the timepoints of interest

n_sig_neurons = length(glm_sig_units);

clear neuron_class
neuron_class.auditory.all = intersect(glm_sig_units,auditory_neuron_idx);
neuron_class.frontal.all = intersect(glm_sig_units,frontal_neuron_idx);

neuron_class.nhp.troy = find(strcmp(spike_log.monkey,'troy'));
neuron_class.nhp.walt = find(strcmp(spike_log.monkey,'walt'));

[sum(strcmp(spike_log.monkey(neuron_class.auditory.all),'troy')), sum(strcmp(spike_log.monkey(neuron_class.auditory.all),'walt'))]
[sum(strcmp(spike_log.monkey(neuron_class.frontal.all),'troy')), sum(strcmp(spike_log.monkey(neuron_class.frontal.all),'walt'))]

%%

for neuron_i = 1:size(spike_log,1)
    anova_id_sig(neuron_i,:) = anova_pvalue{neuron_i}(1,:);
    anova_pos_sig(neuron_i,:) = anova_pvalue{neuron_i}(2,:);
end

figuren('Renderer', 'painters', 'Position', [100 100 850 300]); hold on
id_sub = subplot(1,2,1); hold on
plot(glm_timewin,smooth(nanmean(anova_id_sig(neuron_class.auditory.all,:)),5))
plot(glm_timewin,smooth(nanmean(anova_id_sig(neuron_class.frontal.all,:)),5))
legend({'Auditory','Frontal'}); xlim([-100 600])
title('Identity')
ylabel('P(significant units active)')
pos_sub = subplot(1,2,2); hold on
plot(glm_timewin,smooth(nanmean(anova_pos_sig(neuron_class.auditory.all,:)),5))
plot(glm_timewin,smooth(nanmean(anova_pos_sig(neuron_class.frontal.all,:)),5))
legend({'Auditory','Frontal'}); xlim([-100 600])
title('Position')

linkaxes([id_sub pos_sub],'xy')
subplot(1,2,1); vline([0 563],'k'); vline([413],'k--');
subplot(1,2,2); vline([0 563],'k'); vline([413],'k--');
