function [norm_fr_soundA, norm_fr_soundC, norm_fr_soundG, norm_fr_soundF, norm_fr_soundD, norm_fr_soundAll] =...
    element_extract_normSoundAlign(sdf_soundAlign_data, spike_log)

%% Normalize activity across sounds
% Loop through each neuron in the spike log
for neuron_i = 1:size(spike_log,1)
    fprintf('Normalizing activity for neuron %i of %i... \n', neuron_i, size(spike_log,1));

    % Identify trials for each sound type (A, C, G, F, D)
    a_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'); % & strcmp(string(sdf_soundAlign_data{neuron_i}(:,5)) ,'nonviol');
    c_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'); % & strcmp(string(sdf_soundAlign_data{neuron_i}(:,5)) ,'nonviol');
    g_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'); % & strcmp(string(sdf_soundAlign_data{neuron_i}(:,5)) ,'nonviol');
    f_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'); % & strcmp(string(sdf_soundAlign_data{neuron_i}(:,5)) ,'nonviol');
    d_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'); % & strcmp(string(sdf_soundAlign_data{neuron_i}(:,5)) ,'nonviol');

    % Initialize temporary variable to store the spike density functions for all trials
    clear sound_sdf_temp
    sound_sdf_temp = cell2mat(sdf_soundAlign_data{neuron_i}(~strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1));

    % Calculate the mean and standard deviation of the baseline firing rate from the sound-aligned SDF
    % (using a 100 ms window before sound onset)
    sound_baseline_fr_mu = nanmean(nanmean(sound_sdf_temp(:,200+[-100:0])));
    sound_baseline_fr_std = nanstd(nanmean(sound_sdf_temp(:,200+[-100:0])));

    % Get sound specific SDF
    clear tempSDF_*
    tempSDF_A = cell2mat(sdf_soundAlign_data{neuron_i}(a_trials,1));
    tempSDF_C = cell2mat(sdf_soundAlign_data{neuron_i}(c_trials,1));
    tempSDF_G = cell2mat(sdf_soundAlign_data{neuron_i}(g_trials,1));
    tempSDF_F = cell2mat(sdf_soundAlign_data{neuron_i}(f_trials,1));
    tempSDF_D = cell2mat(sdf_soundAlign_data{neuron_i}(d_trials,1));

    % Normalize firing rates for each sound type by subtracting the mean baseline firing rate
    % and dividing by the standard deviation
    % norm_fr_soundA(neuron_i,:) = smooth((nanmean(tempSDF_A) - sound_baseline_fr_mu) ./ sound_baseline_fr_std,50);
    % norm_fr_soundC(neuron_i,:) = smooth((nanmean(tempSDF_C) - sound_baseline_fr_mu) ./ sound_baseline_fr_std,50);
    % norm_fr_soundG(neuron_i,:) = smooth((nanmean(tempSDF_G) - sound_baseline_fr_mu) ./ sound_baseline_fr_std,50);
    % norm_fr_soundF(neuron_i,:) = smooth((nanmean(tempSDF_F) - sound_baseline_fr_mu) ./ sound_baseline_fr_std,50);
    % norm_fr_soundD(neuron_i,:) = smooth((nanmean(tempSDF_D) - sound_baseline_fr_mu) ./ sound_baseline_fr_std,50);
    norm_fr_soundAll(neuron_i,:) = smooth((nanmean(cell2mat(sdf_soundAlign_data{neuron_i}...
        (find(sum([a_trials, c_trials, g_trials, f_trials, d_trials],2) == 1),1))) - sound_baseline_fr_mu) ./ sound_baseline_fr_std,50);

    norm_fr_soundA(neuron_i,:) = smooth((nanmean(tempSDF_A) - nanmean(nanmean(tempSDF_A(:,200+[-100:0])))) ./ nanstd(nanmean(tempSDF_A(:,200+[-100:0]))),50);
    norm_fr_soundC(neuron_i,:) = smooth((nanmean(tempSDF_C) - nanmean(nanmean(tempSDF_C(:,200+[-100:0])))) ./ nanstd(nanmean(tempSDF_C(:,200+[-100:0]))),50);
    norm_fr_soundG(neuron_i,:) = smooth((nanmean(tempSDF_G) - nanmean(nanmean(tempSDF_G(:,200+[-100:0])))) ./ nanstd(nanmean(tempSDF_G(:,200+[-100:0]))),50);
    norm_fr_soundF(neuron_i,:) = smooth((nanmean(tempSDF_F) - nanmean(nanmean(tempSDF_F(:,200+[-100:0])))) ./ nanstd(nanmean(tempSDF_F(:,200+[-100:0]))),50);
    norm_fr_soundD(neuron_i,:) = smooth((nanmean(tempSDF_D) - nanmean(nanmean(tempSDF_D(:,200+[-100:0])))) ./ nanstd(nanmean(tempSDF_D(:,200+[-100:0]))),50);

end