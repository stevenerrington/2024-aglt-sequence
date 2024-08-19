
%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog, stimulusLog, spike_log] = import_exp_map();
ephysLog = clean_exp_map(ephysLog);
spike_log = clean_spike_map(spike_log);

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
session_i = 34;
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

load(fullfile(dirs.mat_data,outfile_name))

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

clear tdt_data 
tdtFun = @TDTbin2mat;
tdt_data = tdtFun(fullfile(dirs.raw_data,exp_filename),'STORE',{'EyeX','EyeY'});

eye_x = resample(tdt_data.streams.EyeX.data,1000,round(tdt_data.streams.EyeX.fs));
eye_y = resample(tdt_data.streams.EyeY.data,1000,round(tdt_data.streams.EyeY.fs));


ops.timewin = -1000:5000;
timeWin = ops.timewin;

alignTimes = event_table.stimulusOnset_ms;

eye_aligned = nan(length(alignTimes),range(timeWin)+1);

for ii = 1:length(alignTimes)
    try
        eye_aligned(ii,:) = eye_y(alignTimes(ii)+timeWin(1):alignTimes(ii)+timeWin(end));

    end
end


figuren;
plot(ops.timewin, eye_aligned(60,:))

