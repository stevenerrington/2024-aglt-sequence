
%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'troy'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2022-02-09-Opto'; % Experimental raw data
task = 'opto-agl_t'; % Experiment type [agl, opto]
session_n = '0002'; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name



%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

spike_status = 'phy';

switch spike_status
    
    % Convert .ncs files into a binary file for use in Kilosort --------/
    case 'bin' 
        for ch_n = 1:n_channels
            clear filepart_name NCSpath spk_ncs_in

            filepart_name = ['CSC' int2str(ch_n) '_' session_n];
            NCSpath = [fullfile(dirs.raw_data,exp_filename,filepart_name) '.ncs'];

            spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(dirs.raw_data,exp_filename));
            spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
        end

        % Create a binary file and export the restructure broadband data
        clear bin_out_file
        bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
        fwrite(bin_out_file,spk_ncs_out,'int16');
        fclose(bin_out_file);

        % Run kilosort
        mkdir(fullfile(dirs.kilosort,outfile_name));
        % - Run in python: to be integrated here.

    % Import data phy-curated Kilosort data ---------------------------/   
    case 'phy'  
        % Import phy curated data
        ops = struct();
        ops.rootZ = fullfile(dirs.kilosort,outfile_name);
        ops.bin_file = [dirs.bin_data outfile_name '.dat'];
        ops.nCh = n_channels;
        ops.fs = 32000;

        [spikes] = phy2mat(ops);
        [spk_info] = phyinfo2mat(ops);
end

% Local field potential data -------------------------------------------------------
filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1_' session_n '.ncs'],32);
lfp = ft_read_neuralynx_interp(filelabels_lfp);
lfp = lfp.trial{1};

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n; 
clear event_table opto_event
event_table = get_event_table(ops);

% Event extraction
clear opto_event aligntime
ops.event_port = 2;
opto_event = get_opto_trials(event_table, ops);

aligntime = opto_event.stimOn;


% Figure ----------------------------------------------------------------
figuren;
histogram(spikes.time.DSP09a,100)
vline(aligntime,'r')


median(spikes.time.DSP09a)
median(aligntime)

% Get aligned neural activity ------------------------------------------
ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);
lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);



