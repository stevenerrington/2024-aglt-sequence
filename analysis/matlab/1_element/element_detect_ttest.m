function [ttest_encode_flag, ttest_encode_dir] = element_detect_ttest(sdf_soundAlign_data, spike_log)

parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1))
    sound_list = {'A','C','D','F','G'};
    for sound_i = 1:5
        sdf_in = []; sdf_in = cell2mat(sdf_soundAlign_data{neuron_i}(:,1))
        baseline_in = []; baseline_in = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),sound_list{sound_i}),201+[-200:0]),2);
        activity_in = []; activity_in = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),sound_list{sound_i}),200+[0:400]),2);

        ttest_encode_flag(neuron_i,sound_i) = ttest(baseline_in, activity_in);
        ttest_encode_dir(neuron_i,sound_i) = computeCohen_d(baseline_in, activity_in);
    end
end

end