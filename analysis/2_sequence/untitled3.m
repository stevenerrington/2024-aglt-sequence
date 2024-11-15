
seq_rep_n_list = zeros(1,16);

neuron_i = 1;
session_i = find(strcmp(ephysLog.session, spike_log.session{neuron_i}));
monkey =  spike_log.monkey{neuron_i};
area = spike_log.area{neuron_i};

load(fullfile(dirs.mat_data,spike_log.session{neuron_i}),'event_table')

trial = 1;

seq_idx = event_table.cond_value(trial)


seq_pos = 1

if seq_pos == 1
    seq_rep_n_list(seq_idx) = seq_rep_n_list(seq_idx) +1;
    seq_rep_n = seq_rep_n_list(seq_idx);
end


exp_elem = 'A'
confidence = 0.5
act_elem = 'X'
surprise = 'X'
firing_rate = 10.4;

