
%%
%{
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab opto script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%}

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
optoLog = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'opto'));

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
for session_i = 1:size(optoLog,1)
    monkey = optoLog.monkey{session_i}; % Monkey name [troy, chief]

    % Experimental parameters -------------------------------------------
    n_channels = 32; % Number of channels recorded in session

    % Key setup variables
    exp_filename = optoLog.data_folder{session_i}; % Experimental raw data
    task = optoLog.task{session_i}; % Experiment type [agl, opto]
    session_n = optoLog.file_n{session_i}; % Experimental file tag

    % Define experimental/data directories -------------------------------
    outfile_name = optoLog.session{session_i}; % Processed file name

    dirs.raw_data = optoLog.data_dir{session_i};

    try
        % Local field potential data -------------------------------------------------------
        filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n '.ncs'],32);
        lfp = ft_read_neuralynx_interp(filelabels_lfp);
        lfp = lfp.trial{1};

        % Behavioral data -------------------------------------------------------
        % Read in events
        ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n;
        clear event_table_raw event_table
        event_table_raw = get_event_table(ops);
        opto_event = get_opto_trials(event_table_raw);

        aligntime = opto_event.laserOnset_ms;

        ops.timewin = [-1000:5000];
        ops.freq = [1 60];
        ops.ch_extract = [1:32];
        lfp = patch_fault_ch(lfp,23);
        [~, lfp_array] = get_lfp_aligned(lfp,aligntime,ops);
        trial_average_lfp = nanmean(lfp_array,3);

        figuren('Renderer', 'painters', 'Position', [100 100 400 1200]);

        subplot(1,1,1); hold on
        for ch_i = 1:n_channels
            if ismember(ch_i,[1:16])
                color_line = [204 0 102]./255;
            else
                color_line = [0 153 153]./255;
            end

            plot(ops.timewin,  trial_average_lfp(ch_i,:)+10*(ch_i-1),'color',color_line)
        end

        % vline(test,'k-')
        % vline(test2,'r-')
        set(gca,'ydir', 'reverse')
        ylim([-15 (length(1:n_channels)*10)]); yticks([10*([1:32]-1)]); yticklabels(num2cell([1:32]))
        xlim([-250 1000])
        title([outfile_name '-' optoLog.laser_color{session_i} '-'  optoLog.laser_freq{session_i} '-' optoLog.laser_depth{session_i} '-' optoLog.laser_area{session_i}])

        fig = gcf;

        set(fig,'renderer','painters','Units','Inches');
        pos = get(fig,'Position');
        set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(fig,['C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data-extraction\doc\opto_lfp_summary\' outfile_name '.pdf'],'-r400','-bestfit','-dpdf')
        close all
    end
end