function kikuchi_tdt_extractFunction(ephysLog, session_i, dirs)

%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

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

clear tdt_data 
tdtFun = @TDTbin2mat;
tdt_data = tdtFun(fullfile(dirs.raw_data,exp_filename));

% Spike data -------------------------------------------------------
% Preprocess: adjust for variable channel duration between probes
[c, shortest_chan] = min([length(tdt_data.streams.Raws.data) length(tdt_data.streams.Raw1.data)]);

% --! check earlier onsets between channels too
[a, earliest_chan] = min([tdt_data.streams.Raws.startTime tdt_data.streams.Raw1.startTime]);
onset_adj = abs(diff([tdt_data.streams.Raws.startTime tdt_data.streams.Raw1.startTime]));

switch shortest_chan
    case 1
        tdt_data.streams.Raws.data = tdt_data.streams.Raws.data(:,1:length(tdt_data.streams.Raws.data));
        tdt_data.streams.Raw1.data = tdt_data.streams.Raw1.data(:,1:length(tdt_data.streams.Raws.data));
    case 2
        tdt_data.streams.Raws.data = tdt_data.streams.Raws.data(:,1:length(tdt_data.streams.Raw1.data));
        tdt_data.streams.Raw1.data = tdt_data.streams.Raw1.data(:,1:length(tdt_data.streams.Raw1.data));
end

clear spk_ncs_out
spk_ncs_out = [tdt_data.streams.Raws.data; tdt_data.streams.Raw1.data]*1e6;

% Re-referencing signal
spk_ncs_out(1:16,:) = spk_ncs_out(1:16,:) - mean(spk_ncs_out(1:16,:)); % Re-reference to electrode mean
spk_ncs_out(17:32,:) = spk_ncs_out(17:32,:) - mean(spk_ncs_out(17:32,:)); % Re-reference to electrode mean

% Resample signal
for channel = 1:n_channels
    spk_ncs_out_rs(channel,:) = resample(spk_ncs_out(channel,:),20000,round(tdt_data.streams.Raws.fs));
end

% Create a binary file and export the restructure broadband data
clear bin_out_file
bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
fwrite(bin_out_file,spk_ncs_out_rs,'int16');
fclose(bin_out_file);

% Run kilosort
mkdir(fullfile(dirs.kilosort,outfile_name));

% - Run in python: to be integrated here.
spikes = 'SPIKE SORTING REQUIRED'; spk_info = 'SPIKE SORTING REQUIRED';


% LFP (resampled to 1kHz) ///////////////////////////////////////////////
% Tidy and pre-process data
% Match length of recording between channels 
[c, shortest_chan] = min([length(tdt_data.streams.LFP1.data) length(tdt_data.streams.LFP2.data)]);
% --! check earlier onsets between channels too
[a, earliest_chan] = min([tdt_data.streams.LFP1.startTime tdt_data.streams.LFP2.startTime]);
onset_adj = abs(diff([tdt_data.streams.LFP1.startTime tdt_data.streams.LFP2.startTime]));

switch shortest_chan
    case 1
        tdt_data.streams.LFP1.data = tdt_data.streams.LFP1.data(:,1:length(tdt_data.streams.LFP1.data));
        tdt_data.streams.LFP2.data = tdt_data.streams.LFP2.data(:,1:length(tdt_data.streams.LFP1.data));
    case 2
        tdt_data.streams.LFP1.data = tdt_data.streams.LFP1.data(:,1:length(tdt_data.streams.LFP2.data));
        tdt_data.streams.LFP2.data = tdt_data.streams.LFP2.data(:,1:length(tdt_data.streams.LFP2.data));
end

% Extract LFP from raw data and resample
clear lfp
for channel = 1:n_channels
    if ismember(channel,1:16)
        lfp(channel,:) = resample(tdt_data.streams.LFP1.data(channel,:),1000,round(tdt_data.streams.LFP1.fs))*1e6;
    else
        lfp(channel,:) = resample(tdt_data.streams.LFP2.data(channel-16,:),1000,round(tdt_data.streams.LFP2.fs))*1e6;
    end
end

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename;
event_table = get_agl_t_trials_tdt(tdt_data.epocs.St1_,ops);

%% Save extracted data

save(fullfile(dirs.mat_data,[outfile_name '.mat']),'event_table','spikes','spk_info','lfp','-v7.3')
fprintf('Extracted data successfully saved to %s    \n', fullfile(dirs.mat_data,[outfile_name '.mat']))
fprintf(' - Events  ✓   \n')
fprintf(' - Spikes  ✓   \n')
fprintf(' - LFP     ✓   \n')



