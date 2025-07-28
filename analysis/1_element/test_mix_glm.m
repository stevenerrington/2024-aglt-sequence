for area_i = 1:2
    switch area_i
        case 1
            area = 'aud';
        case 2
            area = 'frontal';
    end
    for clu_i = 1:6
        neuron_in = []; neuron_in = neuron_class.cluster_idx.(['clu' int2str(clu_i) '_' area]);
        glm_table = table();
        for neuron_i = 1:length(neuron_in)
            neuron_idx = neuron_in(neuron_i);
            valid_trials = [];
            valid_trials = find(strcmp(sdf_soundAlign_data{neuron_idx}(:,5),'nonviol') & cell2mat(sdf_soundAlign_data{neuron_idx}(:,9)) == 1);
            sdf_in = []; average_fr = [];
            sdf_in = cell2mat(sdf_soundAlign_data{neuron_idx}(:,1));
            average_fr = nanmean(sdf_in(valid_trials,200+[0:413]),2);
            identity = sdf_soundAlign_data{neuron_idx}(valid_trials,3);
            position = sdf_soundAlign_data{neuron_idx}(valid_trials,4);
            neuron_label = repmat(neuron_idx, length(identity),1);
            glm_table = [glm_table; table(average_fr, identity, position, neuron_label)];
        end
        lme1 = fitlme(glm_table, 'average_fr ~ position + (1|neuron_label)');
        lme2 = fitlme(glm_table, 'average_fr ~ identity + (1|neuron_label)');
        bic1 = lme1.ModelCriterion.BIC;
        bic2 = lme2.ModelCriterion.BIC;
        bic_diff(clu_i, area_i) = bic2-bic1;
    end
end

% Assume bic_diff is already computed (6x2 matrix)
areas = {'aud', 'frontal'};
clusters = 1:6;

for area_i = 1:2
    fprintf('--- Area: %s ---\n', areas{area_i});
    for clu_i = clusters
        diff_val = bic_diff(clu_i, area_i);
        if diff_val > 0
            fprintf('Cluster %d: POSITION model is better (BIC diff = %.2f)\n', clu_i, diff_val);
        elseif diff_val < 0
            fprintf('Cluster %d: IDENTITY model is better (BIC diff = %.2f)\n', clu_i, diff_val);
        else
            fprintf('Cluster %d: Both models are equal (BIC diff = 0)\n', clu_i);
        end
    end
    fprintf('\n');
end