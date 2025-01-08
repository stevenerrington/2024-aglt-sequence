
ops.glm_avFR_win = 200+[0:400];
brain_areas.auditory = {'R', 'A1', 'RM', 'dSTS'};
brain_areas.frontal = {'44', '45', 'FOP'};


%%

for neuron_i = 1:size(spike_log,1)

    clear session_i monkey area_label area event_table glm_table
    count = 0;
    seq_rep_n_list = zeros(1,16);

    session_i = find(strcmp(ephysLog.session, spike_log.session{neuron_i}));
    monkey =  spike_log.monkey(neuron_i);
    area_label = spike_log.area(neuron_i);

    if any(strcmp(area_label, brain_areas.auditory))
        area = {'auditory'};
    else
        area = {'frontal'};
    end

    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

    load(fullfile(dirs.mat_data,spike_log.session{neuron_i}),'event_table')

    for trial = 1:size(event_table,1)

        seq_idx = event_table.cond_value(trial);

        if seq_idx < 1 | seq_idx > 16 | isnan(seq_idx)
            seq_idx = 1;
            rwd_trial = 0;
        else
            rwd_trial = ~isnan(event_table.rewardOnset_ms(trial));
        end

        for seq_pos = 1:5
            if seq_pos == 1
                seq_rep_n_list(seq_idx) = seq_rep_n_list(seq_idx) +1;
                seq_rep_n = seq_rep_n_list(seq_idx);
            end

            if seq_pos == 1
                prev_elem = {'X'};
            else
                prev_elem = stimulusLog.(['sound_' int2str(seq_pos-1) '_code'])(seq_idx);
            end

            count = count + 1;
            elem_id = stimulusLog.(['sound_' int2str(seq_pos) '_code'])(seq_idx);

            % Get firing rate from sdf_soundAlign_data table
            trial_cell_idx = find(cell2mat(sdf_soundAlign_data{1, neuron_i}(:,8)) == trial);
            ele_cell_idx = trial_cell_idx(find(strcmp(sdf_soundAlign_data{1, neuron_i}(trial_cell_idx,4),['position_' int2str(seq_pos)])));
            baseline_cell_idx = trial_cell_idx(find(strcmp(sdf_soundAlign_data{1, neuron_i}(trial_cell_idx,4),['position_0'])));

            try
                average_win_fr = nanmean(sdf_soundAlign_data{1, neuron_i}{ele_cell_idx, 1}(ops.glm_avFR_win));
                average_win_spkCount = sum(sdf_soundAlign_data{1, neuron_i}{ele_cell_idx, 2} > 0  &...
                    sdf_soundAlign_data{1, neuron_i}{ele_cell_idx, 2} < 400);

                average_win_delta_spkCount = sum(sdf_soundAlign_data{1, neuron_i}{ele_cell_idx, 2} > 0  &...
                    sdf_soundAlign_data{1, neuron_i}{ele_cell_idx, 2} < 400) -...
                    sum(sdf_soundAlign_data{1, neuron_i}{baseline_cell_idx, 2} > 0  &...
                    sdf_soundAlign_data{1, neuron_i}{baseline_cell_idx, 2} < 400);
            catch
                average_win_fr = NaN;
                average_win_spkCount = NaN;
                average_win_delta_spkCount = NaN;
            end


            glm_table(count,:) = table(neuron_i, session_i, monkey, area, trial, seq_idx, seq_rep_n, seq_pos, elem_id, rwd_trial, prev_elem, average_win_fr, average_win_spkCount,average_win_delta_spkCount);

        end
    end

    writetable(glm_table,fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data-extraction\glm_table',['glm_table_neuron' int2str(neuron_i) '.csv']),'WriteRowNames',true)

end


