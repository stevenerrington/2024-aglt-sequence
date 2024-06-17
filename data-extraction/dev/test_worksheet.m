

%% DEVLEOPMENT

 get_laminar_profile([1:16], lfp, agl_t_event)


test_figure_spk
test_figure_lfp



dirs.raw_data = ops.dirs.raw_data;
raw_filename = ops.filename;
session_n = ops.session_n;

[evnt_raw] = ft_read_event_BA(fullfile(dirs.raw_data,raw_filename,['Events_' session_n '.nev']));


aligntime(1), spikes.time.DSP02a(1)




ops.event_code = 1; ops.exp_type = task;
aligntime = get_event_aligntime(ops,event_table);



[lfp_ncs_in.TimeStamp(1), spk_ncs_in.TimeStamp(1) evnt_raw(1).timestamp]


expmat_dir = 'T:\EPHYS\RAWDATA\NHP\Neuralynx\AGL\Troy\data';
exp_mat = load(fullfile(expmat_dir,'AGL_test_Run#1_05_Jul_2021_11_08_07'));


