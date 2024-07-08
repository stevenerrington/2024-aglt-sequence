%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

clear all; clc
dirs.raw_data = 'T:\EPHYS\RAWDATA\NHP\TDT\WALT\ephys\Tank\';
clear tdt_data 
tdtFun = @TDTbin2mat;
tdt_data = tdtFun(fullfile(dirs.raw_data,'Walter-201130-132046'));

audio_data = tdt_data.streams.Audo.data;

audio_data_downsample = resample(audio_data,1000,round(tdt_data.streams.Audo.fs));

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename;
event_table = get_agl_t_trials_tdt(tdt_data.epocs.St1_,ops);


ops.timewin = [-1000:5000];
audio_aligned = get_lfp_aligned(audio_data_downsample,event_table.stimulusOnset_ms,ops);


figuren; hold on
for sound_i = 1:6
    trials_in = [];
    trials_in = find(event_table.cond_value == sound_i);
    subplot(6,3,sound_i)
    plot(ops.timewin,nanmean(audio_aligned.lfp_1(trials_in,:)))
    vline(0,'k')
end


for sound_i = 7:12
    trials_in = [];
    trials_in = find(event_table.cond_value == sound_i);
    subplot(6,3,sound_i)
    plot(ops.timewin,nanmean(audio_aligned.lfp_1(trials_in,:)))
    vline(0,'k')
end


for sound_i = 13:16
    trials_in = [];
    trials_in = find(event_table.cond_value == sound_i);
    subplot(6,3,sound_i)
    plot(ops.timewin,nanmean(audio_aligned.lfp_1(trials_in,:)))
    vline(0,'k')
end