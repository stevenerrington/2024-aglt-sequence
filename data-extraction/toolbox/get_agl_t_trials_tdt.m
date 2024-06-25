function agl_t_event = get_agl_t_trials_tdt(eventcodes,ops)
clear event_code

event_code.trialStart = 1;
event_code.stimulusOnset = 2;
event_code.rewardOnset = 6;
event_code.rewardOffset = 4;

% Non-violation sequence: 1 to 8
% Violation sequence: 9 to 16
event_code.nonviolCondition = [1:8];
event_code.violCondition = [9:16];


code = eventcodes.data;

trial_start_codes = find(code == 1);


% Trial information
clear trial_n cond_value cond_label
for i = 1:length(trial_start_codes)
    
        % Get trial number
        trial_n(i,1) = i;

        try
            % Condition code
            cond_value(i,1) = code(trial_start_codes(i)+2)-9;
        catch
            cond_value(i,1) = NaN;
        end

            % Condition label
            if ismember(cond_value(i,1), [1:8])
                cond_label{i,1} = 'nonviol';
            elseif ismember(cond_value(i,1),[9:16])
                cond_label{i,1} = 'viol';
            else
                cond_label{i,1} = 'error';
            end


end

    event_labels = fieldnames(event_code);
    event_labels = event_labels(1:4);

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
                    eventcodes.onset...
                    (trial_start_codes(i)+find(trial_codeblock{i,1}...
                    == event_code.(event_labels{event_i}))-1)*1000;
            catch
                event_time_ms.(event_labels{event_i})(i,1) = NaN;
            end

        end
    
end



agl_t_event = table(trial_n, cond_value, cond_label,...
    trial_codeblock,...
    event_time_ms.trialStart,event_time_ms.stimulusOnset,...
    event_time_ms.rewardOnset,event_time_ms.rewardOffset,...
    'VariableNames',{'trial_n','cond_value','cond_label',...
    'trial_codeblock','trialStart_ms','stimulusOnset_ms',...
    'rewardOnset_ms','rewardOffset_ms'});

