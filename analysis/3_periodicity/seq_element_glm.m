function [glm_sequence_out] = seq_element_glm(sound_info_in, transitional_probability)

%% Extract: get relevant data for GLM table

n_neurons = size(sound_info_in,2);
summary_lm_table = table;

for neuron_i = 1:n_neurons

    reg_tbl = table;

    sound_info_in_sdf = cell2mat(sound_info_in{neuron_i}(:,1));

    reg_tbl.sound = string(sound_info_in{neuron_i}(:,3));
    reg_tbl.order_pos = string(sound_info_in{neuron_i}(:,4));
    reg_tbl.condition = string(sound_info_in{neuron_i}(:,5));
    reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';
    reg_tbl.neuron_i =  repmat(neuron_i,length(sound_info_in{neuron_i}(:,3)),1);


    for trial_i = 1:size(reg_tbl)

        if trial_i == 1 | strcmp(reg_tbl.sound(trial_i),'Baseline')
            reg_tbl.backward_prob(trial_i) = NaN;
            reg_tbl.forward_prob(trial_i) = NaN;
        else
            curr_trial_ele = find(strcmp(transitional_probability.elements, reg_tbl.sound(trial_i)));
            prev_trial_ele = find(strcmp(transitional_probability.elements, reg_tbl.sound(trial_i-1)));

            try
                reg_tbl.backward_prob(trial_i) = transitional_probability.backward_prob(prev_trial_ele,curr_trial_ele);
                reg_tbl.forward_prob(trial_i) = transitional_probability.forward_prob(prev_trial_ele,curr_trial_ele);
            catch
                reg_tbl.backward_prob(trial_i) = NaN;
                reg_tbl.forward_prob(trial_i) = NaN;
            end

        end

    end

    clear mean_fr std_fr z_sdf window_sdf window_time n_trials n_times baseline_win_idx analysis_win_idx win_fr
    % Setup spike data into GLM
    mean_fr = nanmean(sound_info_in_sdf(:));
    std_fr = nanstd(sound_info_in_sdf(:));

    for trial_i = 1:size(sound_info_in_sdf,1)
        z_sdf(trial_i,:) = (sound_info_in_sdf(trial_i,:)-mean_fr)./std_fr;
    end

    window_size = 100;
    window_shift = 10;
    [window_sdf, window_time] = movaverage_sdf(z_sdf, window_size, window_shift);

    [n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

    baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
    analysis_win_idx = find(window_time >= 0 & window_time <= 400); % Find the relevant indicies for the timepoints of interest
    win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

    reg_tbl.firing_rate = win_fr; % Add the firing rate over this whole window to the GLM table

    % Summary metrics
    a1 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'A') & strcmp(reg_tbl.condition,'nonviol')));
    a2 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol')));
    a3 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'D') & strcmp(reg_tbl.condition,'nonviol')));
    a4 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'F') & strcmp(reg_tbl.condition,'nonviol')));
    a5 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'G') & strcmp(reg_tbl.condition,'nonviol')));

    b1 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.order_pos,'position_2')));
    b2 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.order_pos,'position_3')));
    b3 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.order_pos,'position_4')));
    b4 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & strcmp(reg_tbl.order_pos,'position_5')));

    c1 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & reg_tbl.backward_prob == 28.6));
    c2 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & reg_tbl.backward_prob == 42.9));

    d1 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & reg_tbl.forward_prob == 50));
    d2 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & reg_tbl.forward_prob == 85.7));
    d3 = nanmean(reg_tbl.firing_rate(strcmp(reg_tbl.sound,'C') & strcmp(reg_tbl.condition,'nonviol') & reg_tbl.forward_prob == 100));

    summary_lm_table(neuron_i,:) = table(neuron_i, a1, a2, a3, a4, a5, b1, b2, b3, b4, c1, c2, d1, d2, d3, 'VariableNames',...
        {'neuron_i','sound_a','sound_c','sound_d','sound_f','sound_g','pos_2','pos_3','pos_4','pos_5','bp_low','bp_high','fp_low','fp_mid','fp_high'});

end


%%

% Restructure data for identity model
identity_data = stack(summary_lm_table, {'sound_a', 'sound_c', 'sound_d', 'sound_f', 'sound_g'}, ...
    'NewDataVariableName', 'firing_rate', ...
    'IndexVariableName', 'sound');
identity_data = identity_data(:, {'neuron_i', 'sound', 'firing_rate'});

% Restructure data for position model
position_data = stack(summary_lm_table, {'pos_2', 'pos_3', 'pos_4', 'pos_5'}, ...
    'NewDataVariableName', 'firing_rate', ...
    'IndexVariableName', 'position');
position_data = position_data(:, {'neuron_i', 'position', 'firing_rate'});

% Restructure data for backwards transition model
bp_data = stack(summary_lm_table, {'bp_low', 'bp_high'}, ...
    'NewDataVariableName', 'firing_rate', ...
    'IndexVariableName', 'backwards');
bp_data = bp_data(:, {'neuron_i', 'backwards', 'firing_rate'});

% Restructure data for forward transition model
fp_data = stack(summary_lm_table, {'fp_low', 'fp_mid', 'fp_high'}, ...
    'NewDataVariableName', 'firing_rate', ...
    'IndexVariableName', 'forwards');
fp_data = fp_data(:, {'neuron_i', 'forwards', 'firing_rate'});


%% Run GLM: trial type
clear glm_output

% Identity model
lme_identity = []; 
identity_data.sound = categorical(identity_data.sound);
identity_data.sound = reordercats(identity_data.sound, {'sound_g', 'sound_f', 'sound_d', 'sound_c','sound_a'});
lme_identity = fitlme(identity_data, 'firing_rate ~ sound + (1|neuron_i)');

% Positional model
lme_positional = []; 
position_data.position = categorical(position_data.position);
position_data.position = reordercats(position_data.position, {'pos_2','pos_3','pos_4','pos_5'});
lme_positional = fitlme(position_data, 'firing_rate ~ position + (1|neuron_i)');

% Backward transition model
lme_backward = []; 
lme_backward = fitlme(bp_data, 'firing_rate ~ backwards + (1|neuron_i)');

% Forward transition model
lme_forward = [];
lme_forward = fitlme(fp_data, 'firing_rate ~ forwards + (1|neuron_i)');


%% Output GLM
glm_sequence_out = table({lme_identity}, {lme_positional}, {lme_backward}, {lme_forward},...
    'VariableNames',{'identity_mdl','position_mdl','backward_mdl','forward_mdl'});

