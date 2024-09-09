% Loop through each session in the ephysLog
for session_i = 1:size(ephysLog,1)

    % Load the session data file
    datafile = ephysLog.session{session_i};
    load(fullfile(dirs.mat_data,datafile),'event_table')
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), datafile)

    % Find all neurons recorded in the current session
    neuron_list = spike_log.unitDSP(strcmp(spike_log.session,datafile));

    % Process each neuron in the session
    for neuron_i = 1:length(neuron_list)
        neuron_label = neuron_list{neuron_i}; % Get the neuron label
        load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data\spike',[datafile '_' neuron_label]))
        sound_align_table = table();

        % Process each trial for the current neuron
        for trial_i = 1:size(event_table,1)


            seq_sdf_in = [];  seq_sdf_in = sdf.sequenceOnset(trial_i,:);
            seq_sdf_detrend = [];
            % Linear regression to detrend
            p = polyfit(ops.timewin, seq_sdf_in, 1);
            trend = polyval(p, ops.timewin);
            seq_sdf_detrend(trial_i, :) = seq_sdf_in - trend;
            trial_table = table();

            if ~strcmp(event_table.cond_label(trial_i),'error') % Skip error trials
                for sound_i = 0:5 % Loop through each sound in the trial

                    if sound_i == 0
                        temp = []; temp = sdf.trialStart(trial_i,zero_offset+[-999:1]);

                        trial = trial_i;
                        element_pos = sound_i;
                        element_label = {'Baseline'};
                        element_id = -999;
                        viol_flag = {'NA'};
                        forward_prob = -999;
                        backward_prob = -999;
                        sdf_table = temp(1:1001);
                        raster_table = {raster.sequenceOnset{trial_i}'};

                    else
                        trial = trial_i;
                        element_pos = sound_i;
                        element_label = stimulusLog.(['sound_' int2str(sound_i) '_code'])(event_table.cond_value(trial_i));
                        element_id = find(ismember(transitional_probability.elements,element_label));

                        switch event_table.cond_label{trial_i}
                            case 'nonviol'
                                viol_flag = {'nonviol'};
                            case 'viol'
                                if sound_i == stimulusLog.violation_pos(event_table.cond_value(trial_i))
                                    viol_flag = {'viol'};
                                else
                                    viol_flag = {'nonviol'};
                                end
                        end

                        try
                            prev_element_idx = find(ismember(transitional_probability.elements,...
                                stimulusLog.(['sound_' int2str(sound_i-1) '_code']){event_table.cond_value(trial_i)}));
                        catch
                            prev_element_idx = 1;
                        end

                        backward_prob = transitional_probability.backward_prob(prev_element_idx,element_id);
                        forward_prob = transitional_probability.forward_prob(prev_element_idx,element_id);
                        sdf_table = seq_sdf_detrend(1,zero_offset+ops.sound_sdf_window+sound_onset_ms(sound_i));
                        raster_table = {raster.sequenceOnset{trial_i}'-sound_onset_ms(sound_i)};


                    end
                    trial_table(sound_i+1,:) = table(trial, element_pos, element_label, element_id, viol_flag, forward_prob, backward_prob, sdf_table, raster_table);
                end
            end
            sound_align_table = [sound_align_table;trial_table];
        end

        save(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data\sound_align', ['sound_align_table_' datafile '_' neuron_label]),'sound_align_table','-v7.3')
    end
end

