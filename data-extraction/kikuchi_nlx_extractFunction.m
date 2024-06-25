function kikuchi_nlx_extractFunction(ephysLog, session_i, dirs)
% Check workspace to see if items are clear. If needed re-run workspace
% setup in kikuchi_main.m

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = ephysLog.monkey{session_i}; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = ephysLog.data_folder{session_i}; % Experimental raw data
task = ephysLog.task{session_i}; % Experiment type [agl, opto]
session_n = ephysLog.file_n{session_i}; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = ephysLog.session{session_i}; % Processed file name

dirs.raw_data = ephysLog.data_dir{session_i};

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

% Convert .ncs files into a binary file for use in Kilosort --------/

for ch_n = 1:n_channels
    clear filepart_name NCSpath spk_ncs_in

    filepart_name = ['CSC' int2str(ch_n) session_n];
    NCSpath = [fullfile(dirs.raw_data,exp_filename,filepart_name) '.ncs'];

    spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(dirs.raw_data,exp_filename));
    spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
end

% Re-referencing signal
spk_ncs_out(1:16,:) = spk_ncs_out(1:16,:) - mean(spk_ncs_out(1:16,:)); % Re-reference to electrode mean
spk_ncs_out(17:32,:) = spk_ncs_out(17:32,:) - mean(spk_ncs_out(17:32,:)); % Re-reference to electrode mean

% Create a binary file and export the restructure broadband data
clear bin_out_file
bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
fwrite(bin_out_file,spk_ncs_out,'int16');
fclose(bin_out_file);

% Run kilosort
mkdir(fullfile(dirs.kilosort,outfile_name));
% - Run in python: to be integrated here.
spikes = 'SPIKE SORTING REQUIRED'; spk_info = 'SPIKE SORTING REQUIRED';


% Local field potential data -------------------------------------------------------
filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n '.ncs'],32);
lfp = ft_read_neuralynx_interp(filelabels_lfp);
lfp = lfp.trial{1};

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n;
clear event_table_raw event_table
event_table_raw = get_event_table(ops);
ops.event_port = 2;
event_table = get_agl_t_trials_nlx(event_table_raw, ops);

%% Save extracted data

save(fullfile(dirs.mat_data,[outfile_name '.mat']),'event_table','spikes','spk_info','lfp','-v7.3')
fprintf('Extracted data successfully saved to %s    \n', fullfile(dirs.mat_data,[outfile_name '.mat']))
fprintf(' - Events  ✓   \n')
fprintf(' - Spikes  ✓   \n')
fprintf(' - LFP     ✓   \n')

end
