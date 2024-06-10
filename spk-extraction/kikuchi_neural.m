clear all; warning off

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'troy'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2021-07-05-09-56-06_Opto'; % Experimental raw data
task = 'opto'; % Experiment type [agl, opto]
session_n = '0011'; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name
set_extract_dirs % Set experimental directories (i.e. data, scripts, etc...)


%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

spike_status = 'phy_import';

switch spike_status
    
    % Convert .ncs files into a binary file for use in Kilosort --------/
    case 'bin' 
        for ch_n = 1:n_channels
            clear filepart_name NCSpath spk_ncs_in

            filepart_name = ['CSC' int2str(ch_n) '_' session_n];
            NCSpath = [fullfile(raw_dir,exp_filename,filepart_name) '.ncs'];

            spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(raw_dir,exp_filename));
            spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
        end

        % Create a binary file and export the restructure broadband data
        clear bin_out_file
        bin_out_file = fopen([bin_dir outfile_name '.dat'],'wb');
        fwrite(bin_out_file,spk_ncs_out,'int16');
        fclose(bin_out_file);

        % Run kilosort
        mkdir(fullfile(ks_dir,outfile_name));
        % - Run in python: to be integrated here.

    % Import data phy-curated Kilosort data ---------------------------/   
    case 'phy_import'  
        % Import phy curated data
        ops = struct();
        ops.rootZ = fullfile(ks_dir,outfile_name);
        ops.bin_file = [bin_dir outfile_name '.dat'];
        ops.nCh = n_channels;
        ops.fs = 32000;

        [spikes] = phy2mat(ops);
        [spk_info] = phyinfo2mat(ops);
end

% Local field potential data -------------------------------------------------------
for ch_n = 1:n_channels
    clear filepart_name NCSpath lfp_ncs_in 

    filepart_name = ['LFP' int2str(ch_n) '_' session_n];
    NCSpath = [fullfile(raw_dir,exp_filename,filepart_name) '.ncs'];

    lfp_ncs_in = readncs([filepart_name '.ncs'],fullfile(raw_dir,exp_filename));
    lfp_ncs_out(ch_n,:) = lfp_ncs_in.dat';
end

% Behavioral data -------------------------------------------------------
% Read in events
ops.raw_dir = raw_dir; ops.filename = exp_filename;
ops.session_n = session_n; 
ops.event_code = 1;

aligntime = get_event_aligntime(ops);


%% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Align data to event
% Define parameters
ops.timewin = -1000:5000;
ops.sdf_filter = 'Gauss';

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);
lfp = get_lfp_aligned(lfp_ncs_out,aligntime,ops);



%% DEVLEOPMENT
raw_dir = ops.raw_dir;
raw_filename = ops.filename;
session_n = ops.session_n;

[evnt_raw] = ft_read_event_BA(fullfile(raw_dir,raw_filename,['Events_' session_n '.nev']));


[lfp_ncs_in.TimeStamp(1), spk_ncs_in.TimeStamp(1) evnt_raw(1).timestamp]


test_figure
expmat_dir = 'T:\EPHYS\RAWDATA\NHP\Neuralynx\AGL\Troy\data';
exp_mat = load(fullfile(expmat_dir,'AGL_test_Run#1_05_Jul_2021_11_08_07'));


