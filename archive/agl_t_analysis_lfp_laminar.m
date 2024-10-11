
ops.timewin = -1000:5000;
ops.freq = [2 200];

for session_i = 1:size(ephysLog,1)

    datafile = ephysLog.session{session_i};
    load(fullfile(dirs.mat_data,datafile))

    % Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Patch faulty channel
    record_idx = find(strcmp(ephysLog.session,datafile),1);

    switch ephysLog.monkey{session_i}
        case 'troy'
            fault_ch_idx = [7 23];
        case 'walt'
            fault_ch_idx = ephysLog.faulty_ch(record_idx);
            if fault_ch_idx == 1
                fault_ch_idx = -999;
            end
    end

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

for session_i = 1:size(ephysLog,1)
    datafile = ephysLog.session{session_i};

    clear laminar_info
    laminar_info = laminar_info_all{session_i};
    get_laminar_plotSummary
    
    % Save figure
    set(fig,'renderer','painters','Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,['C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data-extraction\doc\laminar_summary\' datafile '-laminar.pdf'],'-r400','-bestfit','-dpdf')
    close all
end
