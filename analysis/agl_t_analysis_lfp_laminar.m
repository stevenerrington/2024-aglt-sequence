
ops.timewin = -1000:5000;
ops.freq = [2 200];

for session_i = 23:30

    datafile = ephysLog.session{session_i};
    load(fullfile(dirs.mat_data,datafile))

    % Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Patch faulty channel
    record_idx = find(strcmp(ephysLog.session,datafile),1);
    fault_ch_idx = ephysLog.faulty_ch(record_idx);
    lfp = patch_fault_ch(lfp,fault_ch_idx);

    % Find event information to align data on
    aligntime = event_table.stimulusOnset_ms;

    clear trials*
    trials_input = find(~isnan(aligntime));

    % Run alignment algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

    laminar_info_all{session_i}.auditory = get_laminar_info([1:16], lfp_aligned, event_table);
    laminar_info_all{session_i}.vlpfc = get_laminar_info([17:32], lfp_aligned, event_table);


end

for session_i = 22:27
    clear laminar_info
    laminar_info = laminar_info_all{session_i};
    get_laminar_plotSummary
end
