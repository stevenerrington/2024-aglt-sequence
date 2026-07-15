% Load data
neuron_count = 0; lfp_count = 0; % Initialize counters for neurons and LFPs
zero_offset = 1000;
ops.sound_sdf_window = [-200:800];

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
        neuron_count = neuron_count + 1; % Increment the neuron count
        sound_sdf = {}; % Initialize sound SDF cell array
        count = 0; % Initialize counter for sound SDF

        load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data\spike',[datafile '_' neuron_label]))

        % Process each trial for the current neuron
        for trial_i = 1:size(event_table,1)
            if ~strcmp(event_table.cond_label(trial_i),'error') % Skip error trials
                for sound_i = 0:5 % Loop through each sound in the trial
                    count = count + 1;

                    if sound_i == 0
                        temp = []; temp = repmat(sdf.trialStart(trial_i,zero_offset+[-550:-50]),2,1);
                        sound_sdf{count,1} = temp(1:1001);
                        sound_sdf{count,2} = raster.sequenceOnset{trial_i};
                        sound_sdf{count,3} = 'Baseline';
                        sound_sdf{count,4} = ['position_' int2str(sound_i)];
                        sound_sdf{count,5} = 'Baseline';
                        sound_sdf{count,6} = 'NA';
                        sound_sdf{count,7} = neuron_i;
                        sound_sdf{count,8} = trial_i;
                        sound_sdf{count,9} = ~isnan(event_table.rewardOnset_ms(trial_i));

                    else
                        % Extract and store relevant SDF information for each sound

                        clear sdf_trial
                        % Linear regression to detrend
                        sdf_trial = sdf.sequenceOnset(trial_i,zero_offset+ops.sound_sdf_window+sound_onset_ms(sound_i));
                        % p = polyfit(1:length(ops.sound_sdf_window), sdf.sequenceOnset(trial_i,zero_offset+ops.sound_sdf_window+sound_onset_ms(sound_i)), 1);
                        % trend = polyval(p, 1:length(ops.sound_sdf_window));
                        % sdf_detrended = sdf_trial - trend;


                        sound_sdf{count,1} = sdf_trial;
                        sound_sdf{count,2} = raster.sequenceOnset{trial_i}-sound_onset_ms(sound_i);
                        sound_sdf{count,3} = stimulusLog.(['sound_' int2str(sound_i) '_code']){event_table.cond_value(trial_i)};
                        sound_sdf{count,4} = ['position_' int2str(sound_i)];

                        switch event_table.cond_label{trial_i}
                            case 'nonviol'
                                sound_sdf{count,5} = 'nonviol';
                            case 'viol'
                                if sound_i == stimulusLog.violation_pos(event_table.cond_value(trial_i))
                                    sound_sdf{count,5} = 'viol';
                                else
                                    sound_sdf{count,5} = 'nonviol';
                                end
                        end

                        if sound_i == 1
                            sound_sdf{count,6} = 'NA';
                        else
                            sound_sdf{count,6} = stimulusLog.(['sound_' int2str(sound_i-1) '_code']){event_table.cond_value(trial_i)};
                        end
                        sound_sdf{count,7} = neuron_i;
                        sound_sdf{count,8} = trial_i;
                        sound_sdf{count,9} = ~isnan(event_table.rewardOnset_ms(trial_i));

                    end
                end
            end
        end

        % Store the sound aligned SDF for the current neuron
        neuron_sdfsound_out{neuron_count} = sound_sdf;
        % Calculate and store the baseline activity for the current neuron
        neuron_baseline(neuron_count,:) = nanmean(sdf.trialStart(:,zero_offset+ops.bl_win));
    end
end

sdf_soundAlign_data = neuron_sdfsound_out;
save(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data\sound_align','sdf_soundAlign_data'),'sdf_soundAlign_data','-v7.3')