clear all; clc

% Set directory and filenames
noise_dir = 'T:\projects\Steven\2024-08-07_10-53-37-AudLatencyTest';
filepart_name = ['CSC' int2str(32)];


% Note: I don't know why, but the header in the noise latency test data is
% slightly off, which then causes an error during the conversion. To fix
% it, I added the following to line 191 of readncs.m
%
%               flds{1} = flds{1}(end-6:end);
%
% This then allowed the code to run as normal


% Load signal in
spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(noise_dir));
spk_ncs_out(1,:) = spk_ncs_in.dat';

% Downsample signal so that it's easy to work with
original_fs = 32000;
target_fs = 1000;
[p, q] = rat(target_fs / original_fs, 1e-6);
audio_data_downsample = resample(spk_ncs_out, p, q);

% Get event codes & convert time into seconds/ms
hdr = ft_read_header(get_ncs_filelabel([noise_dir '/'], ['CSC1.ncs'],1));
[event] = ft_read_event_BA(fullfile(noise_dir,['Events.nev']));

for i=1:length(event)
  event(i).sample = (event(i).timestamp-double(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1;
end

for i = 1:length(event)
    event(i).timestampSeconds = (event(i).timestamp - double(hdr.FirstTimeStamp)) *1e-6;
end

% Output event codes as a table so that it's easy (for me) to work with.
event = struct2table(event);
event_table = event;
n_events = size(event_table,1);

for event_i = 1:n_events
    event_table.timestamp_ms(event_i) = (event_table.timestampSeconds(event_i))*1000;
end

% Find unique event codes to try and determine the structure of the code
event_codes = unique(event_table.value);

for event_i = 1:length(event_codes)
    n_events(event_i,1) = sum(event_table.value == event_codes(event_i));
    n_events(event_i,2) = event_codes(event_i);
end

ops.dirs.raw_data = 'T:\projects\Steven\'; ops.filename = '2024-08-07_10-53-37-AudLatencyTest'; ops.session_n = '';
clear event_table_raw event_table
event_table_raw = get_event_table(ops);
ops.event_port = 2;
event_table = get_agl_t_trials_nlx(event_table_raw, ops);

% Once we know which code to align on, get the time in ms.
aligntime = event_table.stimulusOnset_ms;

% Then use these times to align the audio signal on
timewin = [-100:4000];
signal_aligned = nan(length(aligntime),range(timewin)+1);

for ii = 1:length(aligntime)
    try
        signal_aligned(ii,:) = audio_data_downsample(aligntime(ii)+timewin(1):aligntime(ii)+timewin(end));
        signal_aligned_norm(ii,:) = signal_aligned(ii,:)./max(abs(signal_aligned(ii,:)));
        signal_aligned_norm(ii,:) = signal_aligned_norm(ii,:)-mean(signal_aligned_norm(ii,:));
    end
end

% Create a plot to check it looks reasonable.
figuren;
plot(signal_aligned_norm(86,:))
hline(mean(signal_aligned_norm(86,:)))

for trial_i = 1:size(event_table,1)
    try
        onset_latency(trial_i,1) = find(abs(signal_aligned_norm(trial_i,:)) > 0.2,1,'first') - 100;
    catch
        onset_latency(trial_i,1) = NaN;
    end
end

figuren;
histogram(onset_latency,16)
vline(mean(onset_latency),'r-')
vline(median(onset_latency),'r--')