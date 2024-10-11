function [roc_neuron_out_ROC, roc_neuron_out_p, sdf_viol_out, sdf_nonviol_out] = violation_detect_ROC(spike_log, dirs)
 
parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron
    sdf_in = load(fullfile('C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data\spike', [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
    event_table_in = load(fullfile(dirs.mat_data,[spike_log.session{neuron_i} '.mat']),'event_table');

    viol_sdf = []; nonviol_sdf = [];
    viol_sdf = sdf_in.sdf.violation(strcmp(event_table_in.event_table.cond_label,'viol'),:);
    nonviol_sdf = sdf_in.sdf.violation(strcmp(event_table_in.event_table.cond_label,'nonviol'),:);

    analysis_win = [0:400];

    [roc_out] = roc_curve(nanmean(nonviol_sdf(:,1000+analysis_win),1), nanmean(viol_sdf(:,1000+analysis_win),1),'dispt',0);
    [p, ~, ~] = permutationTest(nanmean(nonviol_sdf(:,1000+analysis_win),1), nanmean(viol_sdf(:,1000+analysis_win),1), 100);

    roc_neuron_out_ROC(neuron_i, 1) = roc_out.param.AROC;
    roc_neuron_out_p(neuron_i, 1) = p < 0.01;


    baseline_fr_mu = mean(nanmean(sdf_in.sdf.violation(:,1000+[-100:0])));
    baseline_fr_std = std(nanmean(sdf_in.sdf.violation(:,1000+[-100:0])));
    
    sdf_viol_out(neuron_i,:) = smooth((nanmean(viol_sdf) - baseline_fr_mu) ./ baseline_fr_std,50);
    sdf_nonviol_out(neuron_i,:) = smooth((nanmean(nonviol_sdf) - baseline_fr_mu) ./ baseline_fr_std,50);
end
