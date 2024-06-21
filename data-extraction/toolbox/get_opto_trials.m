function opto_event = get_opto_trials(event_table, ops)

clear event_code

event_code.trialStart = 1;
event_code.stimOn = 9;
event_code.reward = 2;

event_table_code = event_table(event_table.port == ops.event_port, :);
code = dec2bin([event_table_code.ttls]);
code = bin2dec(code(:,2:8));

trial_start_codes = find(code == 1);
event_labels = fieldnames(event_code);

% Trial information
clear trial_n cond_value cond_label
for i = 1:length(trial_start_codes)
    % Get trial number
    trial_n(i,1) = i;
end

% Timestamps
clear trial_codeblock event_time_ms
for i = 1:length(trial_start_codes)

    if i < length(trial_start_codes)
        trial_codeblock{i,1} = code(trial_start_codes(i):trial_start_codes(i+1)-1);
    else
        trial_codeblock{i,1} = code(trial_start_codes(i):end);
    end


    for event_i = 1:length(event_labels)
        try
            event_time_ms.(event_labels{event_i})(i,1) =...
                event_table_code.timestamp_ms...
                (trial_start_codes(i)+find(trial_codeblock{i,1}...
                == event_code.(event_labels{event_i}))-1);
        catch
            event_time_ms.(event_labels{event_i})(i,1) = NaN;
        end

    end

end



opto_event = table(trial_n,trial_codeblock,...
    event_time_ms.trialStart,event_time_ms.stimOn,...
    event_time_ms.reward,...
    'VariableNames',{'trial_n','trial_codeblock','trialStart_ms','stimOn',...
    'Reward'});

