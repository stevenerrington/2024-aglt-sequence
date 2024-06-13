function aligntime = get_event_aligntime(ops, event_table)

raw_dir = ops.raw_dir;
raw_filename = ops.filename;
session_n = ops.session_n;

warning off

if isfield(ops,'event_port')
    event_table = event_table(event_table.port == ops.event_port,:);
end

% Convert ttl pulses from binary to decimal event code

clear marker
switch ops.exp_type
    case 'opto'
        marker = dec2bin([event_table.ttls],8);
        marker = marker(:,8);
        marker = str2num(marker);

    case {'agl_t', 'agl'}
        marker = dec2bin([event_table.ttls]);
        marker = bin2dec(marker(:,2:8));

    case 'vocal'
        marker = dec2bin([event_table.ttls]);
        marker = bin2dec(marker(:,2:8));
end



% Find event code times
idx_marker = find(marker == ops.event_code);
aligntime = event_table.timestamp_ms(idx_marker);


