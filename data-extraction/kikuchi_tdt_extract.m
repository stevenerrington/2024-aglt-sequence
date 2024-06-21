
% Check workspace to see if items are clear. If needed re-run workspace
% setup in kikuchi_main.m

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'walt'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2020-11-02-AGLt'; % Experimental raw data
task = 'agl_t'; % Experiment type [agl, opto]

% Define experimental/data directories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name


%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

tdtFun = @TDTbin2mat;

clear all_data 
all_data = tdtFun(fullfile(dirs.raw_data,exp_filename));


clear lfp_1 lfp_2 
% Local Field Potentials
for channel = 1:size(all_data.streams.LFP1.data,1)
    lfp_1(channel,:) = all_data.streams.LFP1.data(channel,:);
end

for channel = 1:size(all_data.streams.LFP2.data,1)
    lfp_2(channel,:) = all_data.streams.LFP2.data(channel,:);
end



% Spikes

clear spk_ncs_out
channel_i = 0;
for channel = 1:size(all_data.streams.Raw1.data,1)
    channel_i = channel_i + 1;
    spk_ncs_out(channel_i,:) = all_data.streams.Raw1.data(channel,:)*1e6;
end

clear spk_ncs_out
channel_i = 0;
for channel = 1:size(all_data.streams.Raws.data,1)
    channel_i = channel_i + 1;
    spk_ncs_out(channel_i,:) = all_data.streams.Raws.data(channel,:)*1e6;
end


% Create a binary file and export the restructure broadband data
clear bin_out_file
bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
fwrite(bin_out_file,spk_ncs_out,'int16');
fclose(bin_out_file);

% Run kilosort
mkdir(fullfile(dirs.kilosort,outfile_name));


% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n; 
clear event_table events
event_table = get_event_table(ops);


%% Save extracted data

save(fullfile(dirs.mat_data,[outfile_name '.mat']),'event_table','spikes','spk_info','lfp_1','-v7.3')
fprintf('Extracted data successfully saved to %s    \n', fullfile(dirs.mat_data,[outfile_name '.mat']))
fprintf(' - Events  ✓   \n')
fprintf(' - Spikes  ✓   \n')
fprintf(' - LFP     ✓   \n')


