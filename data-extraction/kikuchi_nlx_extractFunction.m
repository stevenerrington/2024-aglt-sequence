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
session_n_agle = ephysLog.file_n_agle{session_i}; % Experimental file tag
session_n_aglt = ephysLog.file_n_aglt{session_i}; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = ephysLog.session{session_i}; % Processed file name

dirs.raw_data = ephysLog.data_dir{session_i};

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

% Convert .ncs files into a binary file for use in Kilosort --------/

for ch_n = 1:n_channels
    clear filepart_name_agle filepart_name_aglt NCSpath_agle NCSpath_aglt spk_ncs_in_agle spk_ncs_in_aglt

    filepart_name_agle = ['CSC' int2str(ch_n) session_n_agle];
    filepart_name_aglt = ['CSC' int2str(ch_n) session_n_aglt];


    NCSpath_agle = [fullfile(dirs.raw_data,exp_filename,filepart_name_agle) '.ncs'];
    NCSpath_aglt = [fullfile(dirs.raw_data,exp_filename,filepart_name_aglt) '.ncs'];

    spk_ncs_in_agle = readncs([filepart_name_agle '.ncs'],fullfile(dirs.raw_data,exp_filename));
    spk_ncs_in_aglt = readncs([filepart_name_aglt '.ncs'],fullfile(dirs.raw_data,exp_filename));

    spk_ncs_out(ch_n,:) = [spk_ncs_in_agle.dat', spk_ncs_in_aglt.dat'];

end

data_times.spk = [size(spk_ncs_in_agle.dat',2), size(spk_ncs_in_aglt.dat',2)];


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
filelabels_lfp_agle = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n_agle '.ncs'],32);
filelabels_lfp_aglt = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' session_n_aglt '.ncs'],32);


lfp_agle = ft_read_neuralynx_interp(filelabels_lfp_agle);
lfp_aglt = ft_read_neuralynx_interp(filelabels_lfp_aglt);

lfp = [lfp_agle.trial{1}, lfp_aglt.trial{1}];

data_times.lfp = [size(lfp_agle.trial{1},2), size(lfp_aglt.trial{1},2)];

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n_agle;
clear event_table_raw_agle event_table_raw_aglt event_table
event_table_raw_agle = get_event_table(ops);
event_table_raw_agle.task = repmat({'AGLe'},size(event_table_raw_agle,1),1);

ops.session_n = session_n_aglt;
event_table_raw_aglt = get_event_table(ops);
event_table_raw_aglt.task = repmat({'AGLt'},size(event_table_raw_aglt,1),1);

event_table_raw_all = [event_table_raw_agle; event_table_raw_aglt];


ops.event_port = 2;
event_table = get_agl_t_trials_nlx(event_table_raw_all, ops);


%% Save extracted data

save(fullfile(dirs.mat_data,[outfile_name '.mat']),'event_table','spikes','spk_info','lfp','data_times','-v7.3')
fprintf('Extracted data successfully saved to %s    \n', fullfile(dirs.mat_data,[outfile_name '.mat']))
fprintf(' - Events  ✓   \n')
fprintf(' - Spikes  ✓   \n')
fprintf(' - LFP     ✓   \n')

end
