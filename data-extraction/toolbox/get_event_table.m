function [event_table] = get_event_table(ops)
warning off

dirs.raw_data = ops.dirs.raw_data;
raw_filename = ops.filename;
session_n = ops.session_n;

% Import events from .nev file
if strcmp(session_n,'')
    hdr = ft_read_header(get_ncs_filelabel(fullfile(dirs.raw_data,[raw_filename '\']), ['CSC1.ncs'],32));
    [event] = ft_read_event_BA(fullfile(dirs.raw_data,raw_filename,['Events.nev']));
else
    hdr = ft_read_header(get_ncs_filelabel(fullfile(dirs.raw_data,[raw_filename '\']), ['CSC1' session_n '.ncs'],32));
    [event] = ft_read_event_BA(fullfile(dirs.raw_data,raw_filename,['Events_' session_n '.nev']));
end

% Convert code timestamps into timestamp seconds
for i=1:length(event)
  event(i).sample = (event(i).timestamp-double(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1;
end

for i = 1:length(event)
    event(i).timestampSeconds = (event(i).timestamp - double(hdr.FirstTimeStamp)) *1e-6;
end

%%///////////////////////////////////////

event = struct2table(event);

% Curate event table to only include codes from relevant ports
%event_table = event(find([event.port] == 2),:);
event_table = event;
n_events = size(event_table,1);

% For each event, calculate the timestamp relative to the first, and then
% convert the timestamp into ms since the first
time_zero = 0;

for event_i = 1:n_events
    event_table.timestamp_ms(event_i) = (event_table.timestampSeconds(event_i)-time_zero)*1000;
end


