function [ttest_encode_flag, ttest_encode_dir] = element_detect_ttest(sdf_soundAlign_data, spike_log)

parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1))

        sdf_in = []; sdf_in = cell2mat(sdf_soundAlign_data{neuron_i}(:,1))
        baseline_in = []; baseline_in = nanmean(sdf_in(~strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),201+[-100:0]),2);
        activity_in = []; activity_in = nanmean(sdf_in(~strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),200+[50:250]),2);

        ttest_encode_flag(neuron_i,1) = ttest(baseline_in, activity_in,'alpha',0.01);
        ttest_encode_dir(neuron_i,1) = computeCohen_d(baseline_in, activity_in);

end

end

