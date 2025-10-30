
%% Run GLM
% Parallel loop to perform GLM analysis for each neuron
if ~exist(fullfile(dirs.root,'data','glm_modulation_data.mat'))
    clear glm_encoding_flag glm_beta glm_sig anova_pvalue window_sdf window_time
    parfor neuron_i = 1:size(spike_log,1)
        fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

        % Perform GLM analysis on sound-aligned SDF data for the current neuron
        [glm_encoding_flag(neuron_i,:), glm_beta{neuron_i}, glm_sig{neuron_i},...
            window_sdf{neuron_i,1}, window_time(neuron_i,:)] = glm_sound_modulation(sdf_soundAlign_data{neuron_i});
        [ttest_encoding_flag(neuron_i,:)] = ttest_sound_modulation(sdf_soundAlign_data{neuron_i});
    end
    save(fullfile(dirs.root,'data','glm_modulation_data.mat'),'glm_encoding_flag','glm_beta','glm_sig','window_sdf','window_time','ttest_encoding_flag','-v7.3')
else
    load(fullfile(dirs.root,'data','glm_modulation_data.mat'));
end

% Identify neurons with significant modulation based on GLM results
glm_sig_units = find(sum(glm_encoding_flag, 2) > 0); % Find neurons with significant encoding in any sound category

% For future use: state the time windows for the GLM
glm_timewin = window_time(1,:); % Define the time window for plotting beta weights
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

clear anova_id_sig glm_neuron_sigxtime_a glm_neuron_sigxtime_b
for neuron_i = 1:size(spike_log,1)
    glm_neuron_sigxtime_a(neuron_i,:) = sum(glm_sig{neuron_i}([1:5],:)) > 0;
    glm_neuron_sigxtime_b(neuron_i,:) = sum(glm_sig{neuron_i}([6:10],:)) > 0;
    %glm_neuron_etaxtime_a(neuron_i,:) = anova_pvalue{neuron_i}(3,:);
    %glm_neuron_etaxtime_b(neuron_i,:) = anova_pvalue{neuron_i}(4,:);

end

figuren('Renderer', 'painters', 'Position', [100 100 400 300]); hold on
id_sub = subplot(1,1,1); hold on
plot(glm_timewin,smooth(nanmean(glm_neuron_sigxtime_a(neuron_class.auditory.all,:)),5)','r-')
plot(glm_timewin,smooth(nanmean(glm_neuron_sigxtime_a(neuron_class.frontal.all,:)),5)','b-')

plot(glm_timewin,smooth(nanmean(glm_neuron_sigxtime_b(neuron_class.auditory.all,:)),5)','r--')
plot(glm_timewin,smooth(nanmean(glm_neuron_sigxtime_b(neuron_class.frontal.all,:)),5)','b--')


legend({'Auditory','Frontal'}); xlim([-100 600])
ylabel('P(significant units active)')
vline([0 563],'k'); vline([413],'k--');


%%
neuron_class.auditory.all = intersect(glm_sig_units,auditory_neuron_idx);
neuron_class.frontal.all = intersect(glm_sig_units,frontal_neuron_idx);

figuren('Renderer', 'painters', 'Position', [680,458,169,420]);
subplot(2,1,1)
donut([length(neuron_class.auditory.all) length(auditory_neuron_idx)-length(neuron_class.auditory.all)])

subplot(2,1,2)
donut([length(neuron_class.frontal.all) length(auditory_neuron_idx)-length(neuron_class.frontal.all)])

%%

metric = mean(normFR_in.norm_fr_soundAll(glm_sig_units, 200:613), 2);
[~, sortIdx] = sort(metric);

colorscale = abs(flipud(cbrewer2('seq', 'RdBu', 100)));

figuren('Renderer', 'painters', 'Position', [680,458,250,420]);
imagesc(-200:800,1:length(sortIdx),normFR_in.norm_fr_soundAll(glm_sig_units(sortIdx),:))
xlim([-100 800]); ylim([1 length(sortIdx)]); clim([-3 3])
vline([0 413],'k-')
vline([563],'k--')
colormap(colorscale)

