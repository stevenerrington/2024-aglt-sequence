function aligntime = get_event_aligntime(ops)

raw_dir = ops.raw_dir;
raw_filename = ops.filename;
session_n = ops.session_n;

warning off


% Import events from .nev file
[evnt_raw] = ft_read_event_BA(fullfile(raw_dir,raw_filename,['Events_' session_n '.nev']));
event = struct2table(evnt_raw);

% Curate event table to only include codes from relevant ports
event_table = find([event.port] == 2);
event_table = event(event_table,:);

n_events = size(event_table,1);

for event_i = 1:n_events
    event_table.timestamp_zero(event_i) = (event_table.timestamp(event_i) - event.timestamp(1));
    event_table.timestamp_ms(event_i) = (event_table.timestamp_zero(event_i)/100000)*1000;
end


% Convert ttl pulses from binary to decimal event code
marker = dec2bin([event_table.ttls]);
marker = bin2dec(marker(:,2:8));

% Find event code times
trlbegin_idx = find(marker == ops.event_code);
aligntime = event_table.timestamp_ms(trlbegin_idx);

