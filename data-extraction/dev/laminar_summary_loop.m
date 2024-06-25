
%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

%% Find task files in raw directory
dirs.raw_data = 'T:\EPHYS\RAWDATA\NHP\Neuralynx\AGL\Troy\RAWDATA';

d = dir(dirs.raw_data);
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
dfolders = struct2table(dfolders);

filename_list = dfolders.name;

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'troy'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

for file_i = 41
    % Key setup variables
    try
        exp_filename = filename_list{file_i}; % Experimental raw data
        task = 'agl_t'; % Experiment type [agl, opto]
        session_n = '0004'; % Experimental file tag

        % Define experimental/data directories -------------------------------
        outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name

        %% Extract neurophysiology data
        % Local field potential data -------------------------------------------------------
        filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1_' session_n '.ncs'],32);
        lfp = ft_read_neuralynx_interp(filelabels_lfp);
        lfp = lfp.trial{1};

        %% Extract behavioral/task data
        % Read in events
        ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n;
        clear event_table events
        event_table = get_event_table(ops);

        %% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find event information to align data on
        ops.event_port = 2;
        events = get_agl_t_trials(event_table, ops);
        aligntime = events.stimulusOnset_ms;

        % Run alignment algorithms
        ops.timewin = -1000:5000;
        lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);

        %% Laminar analyses
        laminar_info.auditory = get_laminar_info([1:16], lfp_aligned, events);
        laminar_info.vlpfc = get_laminar_info([17:22,24:32], lfp_aligned, events);

        get_laminar_plotSummary
    end
end