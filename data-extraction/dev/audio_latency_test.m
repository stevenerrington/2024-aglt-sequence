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

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

dirs.raw_data = 'T:\EPHYS\RAWDATA\NHP\TDT\WALT\ephys\Tank\';
clear tdt_data 
tdtFun = @TDTbin2mat;

for session_i = 1:14
    tdt_data = tdtFun(fullfile(dirs.raw_data,ephysLog_hpc.data_folder{session_i}),'STORE', {'Audo','St1_'});

    audio_data = tdt_data.streams.Audo.data;
    original_fs = tdt_data.streams.Audo.fs;
    target_fs = 1000;


    % Calculate the resampling ratio
    [p, q] = rat(target_fs / original_fs, 1e-6);

    % Resample the signal
    audio_data_downsample = resample(audio_data, p, q);

    % Behavioral data -------------------------------------------------------
    % Read in events
    ops.dirs.raw_data = dirs.raw_data; ops.filename = ephysLog_hpc.data_folder{session_i};
    event_table = get_agl_t_trials_tdt(tdt_data.epocs.St1_,ops);

    ops.timewin = [-1000:5000];
    audio_aligned = get_lfp_aligned(audio_data_downsample,event_table.stimulusOnset_ms,ops);

    onset_latency = [];


    for trial_i = 1:size(event_table,1)
        try
            onset_latency(trial_i,1) = find(abs(audio_aligned.lfp_1(trial_i,:)) > 0.02,1,'first') - 1000;
        catch
            onset_latency(trial_i,1) = NaN;
        end
    end


    session_audio_latency{session_i,1} = onset_latency;

end


% 
% 
% figure;
% plot(onset_latency)
% 
% figure;
% subplot(3,1,1)
% plot(audio_aligned.lfp_1(1,:))
% subplot(3,1,2)
% plot(audio_aligned.lfp_1(90,:))
% subplot(3,1,3)
% plot(audio_aligned.lfp_1(150,:))
% 
% 
% 
% figuren; hold on
% for sound_i = 1:6
%     trials_in = [];
%     trials_in = find(event_table.cond_value == sound_i);
%     subplot(6,3,sound_i)
%     plot(ops.timewin,nanmean(audio_aligned.lfp_1(trials_in,:)))
%     vline(0,'k')
% end
% 
% 
% for sound_i = 7:12
%     trials_in = [];
%     trials_in = find(event_table.cond_value == sound_i);
%     subplot(6,3,sound_i)
%     plot(ops.timewin,nanmean(audio_aligned.lfp_1(trials_in,:)))
%     vline(0,'k')
% end
% 
% 
% for sound_i = 13:16
%     trials_in = [];
%     trials_in = find(event_table.cond_value == sound_i);
%     subplot(6,3,sound_i)
%     plot(ops.timewin,nanmean(audio_aligned.lfp_1(trials_in,:)))
%     vline(0,'k')
% end