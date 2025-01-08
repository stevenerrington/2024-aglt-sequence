function roc_neuron_out = element_detect_ROC_all(sdf_soundAlign_data,spike_log)

for neuron_i = 1:size(spike_log,1)
    fprintf('Running ROC and permutation analysis for neuron %i of %i... \n', neuron_i, size(spike_log,1));

    sdf_in = []; sdf_in = cell2mat(sdf_soundAlign_data{neuron_i}(:,1));

    clear fr*
    fr_baseline = nanmean(sdf_in(~strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,[400:800]),2);


    fr_element = nanmean(sdf_in(~strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[0:400]),2);


    [roc_out] = roc_curve(fr_baseline, fr_element);
    [p, ~, ~] = permutationTest(fr_baseline, fr_element, 100);

    roc_neuron_out(neuron_i, 1) = roc_out.param.AROC;
    roc_neuron_out(neuron_i, 2) = p < 0.05;

end
end
