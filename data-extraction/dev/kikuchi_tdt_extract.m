
%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Check workspace to see if items are clear. If needed re-run workspace
% setup in kikuchi_main.m

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'walt'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2020-11-30-AGLt'; % Experimental raw data
task = 'agl_t'; % Experiment type [agl, opto]

% Define experimental/data directories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

clear tdt_data 
tdtFun = @TDTbin2mat;
tdt_data = tdtFun(fullfile(dirs.raw_data,exp_filename));

% Neural data -----------------------------------------------------------
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
clear lfp_out
for channel = 1:n_channels
    if ismember(channel,1:16)
        lfp_out(channel,:) = resample(tdt_data.streams.LFP1.data(channel,:),1000,round(tdt_data.streams.LFP1.fs))*1e6;
    else
        lfp_out(channel,:) = resample(tdt_data.streams.LFP2.data(channel-16,:),1000,round(tdt_data.streams.LFP2.fs))*1e6;
    end
end

lfp = lfp_out;


% Spiking activity ///////////////////////////////////////////////////////

spike_status = 'phy';

switch spike_status
    case 'bin'
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

        % % Resample signal
        % for channel = 1:n_channels
        %     spk_ncs_out_rs(channel,:) = resample(spk_ncs_out(channel,:),24000,round(tdt_data.streams.Raws.fs));
        % end

        spk_ncs_out_rs = spk_ncs_out;

        % Create a binary file and export the restructure broadband data
        clear bin_out_file
        bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
        fwrite(bin_out_file,spk_ncs_out_rs,'int16');
        fclose(bin_out_file);

        % Run kilosort
        mkdir(fullfile(dirs.kilosort,outfile_name));

    case 'phy'
        % Import phy curated data
        ops = struct();
        ops.rootZ = fullfile(dirs.kilosort,outfile_name);
        ops.bin_file = [dirs.bin_data outfile_name '.dat'];
        ops.nCh = n_channels;
        ops.fs = 24414.0625;

        [spikes] = phy2mat(ops);
        [spk_info] = phyinfo2mat(ops);

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



