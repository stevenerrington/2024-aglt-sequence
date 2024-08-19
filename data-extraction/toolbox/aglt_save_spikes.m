%% Extract sound aligned Spike Density Function (SDF)
% Load data
neuron_count = 0; lfp_count = 0; % Initialize counters for neurons and LFPs

% Loop through each session in the ephysLog
for session_i = 1:size(ephysLog,1)

    % Load the session data file
    datafile = ephysLog.session{session_i};
    load(fullfile(dirs.mat_data,datafile))
    fprintf('Session %i of %i | %s \n', session_i, size(ephysLog,1), datafile)

    % Adjust stimulus onset times based on session-specific audio latency measures
    for trial_i = 1:size(event_table,1)
        event_table.stimulusOnset_ms(trial_i) = event_table.stimulusOnset_ms(trial_i) + session_audio_latency{session_i}(trial_i,1) - 20;
    end

    % Create a violation alignment event for SDF based on condition value
    for trial_i = 1:size(event_table,1)
        switch event_table.cond_value(trial_i)
            case {3, 7, 14} % Conditions with violation time of 1127 ms
                event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 1127;
            case {4, 8, 15} % Conditions with violation time of 2253 ms
                event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 2253;
            case {1, 5, 13} % Conditions with violation time of 1127 ms
                event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 1127;
            case {2, 6, 16} % Conditions with violation time of 2253 ms
                event_table.violation_ms(trial_i) =  event_table.stimulusOnset_ms(trial_i) + 2253;
            case {9, 10, 11, 12} % Conditions without violation
                event_table.violation_ms(trial_i) =  NaN;
        end
    end

    % Find all neurons recorded in the current session
    neuron_list = spike_log.unitDSP(strcmp(spike_log.session,datafile));

    clear sdf_* raster_* aligntime_*
    % Align spikes and generate SDF for different event types
    % - trial start alignment
    aligntime_trialStart = event_table.trialStart_ms;
    [sdf_trialstart, raster_trialstart] = get_spikes_aligned(spikes,aligntime_trialStart,ops);
    % - sequence onset alignment
    aligntime_sequence = event_table.stimulusOnset_ms;
    [sdf_sequence, raster_sequence] = get_spikes_aligned(spikes,aligntime_sequence,ops);
    % - violation onset alignment
    aligntime_viol = event_table.violation_ms;
    [sdf_viol, raster_viol] = get_spikes_aligned(spikes,aligntime_viol,ops);
    % - reward onset alignment
    aligntime_reward = event_table.rewardOnset_ms;
    [sdf_rwd, raster_rwd] = get_spikes_aligned(spikes,aligntime_reward,ops);

    % Process each neuron in the session
    for neuron_i = 1:length(neuron_list)
        neuron_label = neuron_list{neuron_i}; % Get the neuron label
        neuron_count = neuron_count + 1; % Increment the neuron count
        count = 0; % Initialize counter for sound SDF

        clear sdf raster
        sdf.trialStart = sdf_trialstart.(neuron_label);
        sdf.sequenceOnset = sdf_sequence.(neuron_label);
        sdf.violation = sdf_viol.(neuron_label);
        sdf.reward = sdf_rwd.(neuron_label);

        raster.trialStart = raster_trialstart.(neuron_label);
        raster.sequenceOnset = raster_sequence.(neuron_label);
        raster.violation = raster_viol.(neuron_label);
        raster.reward = raster_rwd.(neuron_label);

        save(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-laminar\data\spike',[datafile '_' neuron_label]),'sdf','raster','-v7.3')

    end
end