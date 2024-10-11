%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- AGLt multi-session extraction script  -------------------------------
      S P Errington, 2024

This script was developed with the intention to loop through raw data files
listed in the ephysLog, determine what system is was recorded with, and
then extract the relevant signals and information. This should allow for a
common data format (spikes, spk info, event times, and lfps) to be
extracted for data recorded on both NeuraLynx and TDT.
///////////////////////////////////////////////////////////////////////////
%} 

%% Run initial extraction
% 

% Get session number
for session_i = 1:size(ephysLog,1)
    % Clear workspace to reduce clutter
    clearvars -except C session_i dirs ephysLog

    % Detemine recording system
    exp_system = ephysLog.sys{session_i}; % Experiment type [agl, opto]

    % Provide user update
    h = waitbar(session_i/size(ephysLog,1),{['Processing session ' int2str(session_i) ' of ' int2str(size(ephysLog,1))], ephysLog.session{session_i}});

    % Run extraction functions
    switch exp_system
        case 'plex' % For Neuralynx
            kikuchi_nlx_extractFunction(ephysLog, session_i, dirs)
        case 'tdt' % For TDT
            kikuchi_tdt_extractFunction(ephysLog, session_i, dirs)
    end
    
    close(h)
end

%% Import spike sorted data and add to processed data file
% Once Kilosort has been run

for session_i = 1:size(ephysLog,1)
    outfile_name = ephysLog.session{session_i}; % Processed file name
    exp_system = ephysLog.sys{session_i}; % Experiment type [agl, opto]

    switch exp_system
        case 'plex' % For Neuralynx (32000Hz sampling)
            kikuchi_phy_import(outfile_name, dirs, 32000)
        case 'tdt' % For TDT  (24414.0625Hz sampling)
            kikuchi_phy_import(outfile_name, dirs, 24414.0625)
    end

end


%% Create a spike/unit summary log

% Import and curate experimental log
ephysLog_all = import_exp_map();
ephysLog = clean_exp_map(ephysLog_all);

spike_log = [];

for session_i = 1:size(ephysLog,1)
    outfile_name = ephysLog.session{session_i}; % Processed file name
    load(fullfile(dirs.mat_data,[outfile_name '.mat']),'spk_info');

    % Admin information
    n_units_session = size(spk_info,1);
    n_units_probe1 = sum(ismember(spk_info.site,[1:16]));
    n_units_probe2 = sum(ismember(spk_info.site,[17:32]));
    
    session = repmat(ephysLog.session(session_i),n_units_session,1);
    monkey = repmat(ephysLog.monkey(session_i),n_units_session,1);

    session_allidx = find(strcmp(ephysLog.session{session_i},ephysLog_all.session));
        
    area = [repmat(ephysLog_all.area_label_recA(session_allidx(1)),n_units_probe1,1);...
        repmat(ephysLog_all.area_label_recA(session_allidx(2)),n_units_probe2,1)];

    ml = repmat(ephysLog.nan_x(session_i),n_units_session,1);
    ap = repmat(ephysLog.nan_y(session_i),n_units_session,1);

    admin_table = table(session, monkey, area, ml, ap);
    spike_log = [spike_log; [admin_table, spk_info]];
end

for spike_i = 1:size(spike_log)
    spike_log.high_isi_viol(spike_i,:) = spike_log.ISI_2ms(spike_i) > 2;
    spike_log.low_spk_viol(spike_i,:) = spike_log.nSpikes(spike_i) < 1000;
end

% Save table from import
writetable(spike_log,fullfile(dirs.doc_data,'spike_log.csv'),'WriteRowNames',true)  

% Upload to Google Sheet if required

% Now reference Google sheet
spike_log = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t_spikes'));

%% Extract laminar information

ops.timewin = -1000:5000;
ops.freq = [2 100];

for session_i = 1:size(ephysLog,1)

    datafile = ephysLog.session{session_i};
    data_in = load(fullfile(dirs.mat_data,datafile));
    
    fprintf('Session %i of %i \n', session_i, size(ephysLog,1))
    % Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Patch faulty channel
    record_idx = find(strcmp(ephysLog.session,datafile),1);
    fault_ch_idx = ephysLog.faulty_ch(record_idx);

    switch ephysLog.monkey{session_i}
        case 'troy'
            lfp = patch_fault_ch(data_in.lfp,[7 23]);
        case 'walt'
            lfp = patch_fault_ch(data_in.lfp,[9]);
    end


    % Find event information to align data on
    aligntime = data_in.event_table.stimulusOnset_ms;

    trials_input = find(~isnan(aligntime));

    % Run alignment algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

    laminar_info_all{session_i}.auditory = get_laminar_info([1:16], lfp_aligned, data_in.event_table);
    laminar_info_all{session_i}.vlpfc = get_laminar_info([17:32], lfp_aligned, data_in.event_table);

end

for session_i = 1:size(ephysLog,1)
    datafile = ephysLog.session{session_i};

    clear laminar_info
    laminar_info = laminar_info_all{session_i};
    get_laminar_plotSummary
    fig = gcf;

    set(fig,'renderer','painters','Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,['C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data-extraction\doc\laminar_summary\' datafile '.pdf'],'-r400','-bestfit','-dpdf')
    close all
end
