function roc_neuron_out = element_detect_ROC(sdf_soundAlign_data,spike_log)

for neuron_i = 1:size(spike_log,1)
    fprintf('Running ROC and permutation analysis for neuron %i of %i... \n', neuron_i, size(spike_log,1));

    sdf_in = []; sdf_in = cell2mat(sdf_soundAlign_data{neuron_i}(:,1));

    clear fr*
    fr_baseline = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,[400:800]),2);


    fr_A_bl = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[-100:0]),2);
    fr_C_bl = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[-100:0]),2);
    fr_D_bl = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[-100:0]),2);
    fr_F_bl = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[-100:0]),2);
    fr_G_bl = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[-100:0]),2);

    fr_A = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[0:400]),2);
    fr_C = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[0:400]),2);
    fr_D = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[0:400]),2);
    fr_F = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[0:400]),2);
    fr_G = nanmean(sdf_in(strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) == 1,200+[0:400]),2);

    fr_bl = {fr_A_bl, fr_C_bl, fr_D_bl, fr_F_bl, fr_G_bl};
    fr_all = {fr_A, fr_C, fr_D, fr_F, fr_G};

    for sound_i = 1:5

        [roc_out] = roc_curve(fr_baseline, fr_all{sound_i});
        [p, ~, ~] = permutationTest(fr_baseline, fr_all{sound_i}, 100);

        roc_neuron_out(neuron_i, sound_i, 1) = roc_out.param.AROC;
        roc_neuron_out(neuron_i, sound_i, 2) = p < 0.01;
    end
end

end