%{
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%}

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off
catch_counter = 0;

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog, stimulusLog, spike_log] = import_exp_map();
ephysLog = clean_exp_map(ephysLog);

%% Setup environment
% Admin --------------------------------------------------------------
for session = 1:size(ephysLog,1)
    monkey = ephysLog.monkey{session}; % Monkey name [troy, chief]
    outfile_name = ephysLog.session{session}; % Processed file name
    session_label = ephysLog.file_n{session};

    % Experimental parameters -------------------------------------------
    n_channels = 32; % Number of channels recorded in session

    data_dir = ephysLog.data_dir{session};
    data_file = ephysLog.data_folder{session};

    %% Extract broadband signal
    clear spk_ncs_out

    switch ephysLog.sys{session}
        case 'plex'
            for ch_n = 1:n_channels
                clear filepart_name NCSpath spk_ncs_in

                filepart_name = ['CSC' int2str(ch_n) session_label];
                NCSpath = [fullfile(data_dir,data_file,filepart_name) '.ncs'];

                spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(data_dir,data_file));
                spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
            end
        case 'tdt'
            for ch_n = 1:n_channels
                % INSERT HERE!
            end
    end

    switch ephysLog.sys{session}
        case 'plex'
            ops.fs = 32000;
        case 'tdt'
            ops.fs = 24414.0625;
    end

    %% Run HOSD Spike Detection Algorithm

    % Define algorithm parameters
    params.lowpass = 4000; % lowpass cutoff ( Hz )
    params.highpass = 200; % highpass cutoff ( Hz )
    params.limit_memory = 5e8;
    params.check_features = false;
    params.randseed = 1;
    saveAllRaw = 0; % boolean to keep spikes data (rather than just the clustered data)
    saveFigFile = 0;

    % Setup data for function
    clear data
    data.dat = spk_ncs_out;
    data.fs = ops.fs;

    % - 2024-06-14: Ran this code and it took an age (stopped after 45 minutes
    % to access MATLAB for other tasks). Will run again over the weekend?
    mat_savefile = fullfile('C:\KIKUCHI-LOCAL\data\ephys\hosd',[outfile_name '-HOSD.mat']);

    try
        if exist(mat_savefile) ~= 2
            [spikes,hos,znrm] = HOSD_spike_detection(data,params); % Run extraction
            spike_cluster = sort_spikes(spikes); % Sort spikes
            save(mat_savefile,'spikes','hos','znrm','spike_cluster','-v7.3')
        end
    catch
        catch_counter = catch_counter + 1;
        error_log{catch_counter, 1} = mat_savefile;
    end

end