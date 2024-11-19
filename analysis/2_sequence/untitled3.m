
seq_rep_n_list = zeros(1,16);
ops.glm_avFR_win = 200+[0:400];
transitional_probability = estimate_transitional_probability();

%%
neuron_i = 1;
session_i = find(strcmp(ephysLog.session, spike_log.session{neuron_i}));
monkey =  spike_log.monkey(neuron_i);
area = spike_log.area(neuron_i);

load(fullfile(dirs.mat_data,spike_log.session{neuron_i}),'event_table')

%%
% Initialize transition probability matrix with learned prior probabilities
% Assuming we have a prior transition probability matrix based on compatible sequences
% Example: Replace with actual prior probabilities if available
sounds = {'Start', 'YAG', 'KEM', 'LEK', 'PAV', 'RAP', 'X'};  % Define your sounds

for seq_i = 1:16
    observedSequenceList{seq_i}{1,1} = 'Start';
    for ele_i = 2:6
    observedSequenceList{seq_i}{1,ele_i} =...
        stimulusLog.(['sound_' int2str(ele_i-1) '_label']){seq_i};
    end
    observedSequenceList{seq_i}{1,7} = 'X';
end

for trial_i = 1:size(event_table,1)
    trial_seq = event_table.cond_value(trial_i);
    observedSequence_trial{trial_i,1} = observedSequenceList{trial_seq};
end

%%

trial = 1;

seq_idx = event_table.cond_value(trial);

seq_pos = 1;

if seq_pos == 1
    seq_rep_n_list(seq_idx) = seq_rep_n_list(seq_idx) +1;
    seq_rep_n = seq_rep_n_list(seq_idx);
    exp_elem = {'YAG'};
    confidence = 1;
    act_elem = {'YAG'};
    surprise = 0;
end

seq_pos = 2;
exp_elem = {'A'};



transitional_probability.forward_prob

confidence = 1;
act_elem = observedSequence_trial{trial_i,1}(seq_pos+1);
surprise = 0;







% Get firing rate from sdf_soundAlign_data table
trial_cell_idx = find(cell2mat(sdf_soundAlign_data{1, neuron_i}(:,8)) == trial);
ele_cell_idx = find(strcmp(sdf_soundAlign_data{1, neuron_i}(trial_cell_idx,4),['position_' int2str(seq_pos)]));

average_win_fr = nanmean(sdf_soundAlign_data{1, neuron_i}{ele_cell_idx, 1}(ops.glm_avFR_win));




test = table(neuron_i, session_i, monkey, area, trial, seq_idx, seq_rep_n, seq_pos, exp_elem, confidence, act_elem, surprise, average_win_fr)


