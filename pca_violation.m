% Clear previously stored data variables related to SDFs
clear viol_sdf nonviol_sdf

% Loop over each neuron in spike_log
for neuron_i = 1:size(spike_log,1)

    comp_list = {'G','position_3'; 'F', 'position_5'; 'F', 'position_3'; 'C', 'position_5'};

    % Compute baseline firing rate (mean and std) for z-scoring
    baseline_mu_fr = nanmean(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));
    baseline_std_fr = nanstd(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));

    viol_sdf_loop = []; nonviol_sdf_loop = [];
    next_viol_sdf_loop = []; next_nonviol_sdf_loop = [];

    for comp_i = 1:size(comp_list,1)

        valid_trials = []; nonvalid_trials = [];

        % Identify valid trials: correct and labeled as 'nonviol'
        valid_trials = cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'nonviol') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),comp_list{comp_i,1}) & strcmp(sdf_soundAlign_data{neuron_i}(:,4),comp_list{comp_i,2});
        nonvalid_trials = cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'viol') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),comp_list{comp_i,1}) & strcmp(sdf_soundAlign_data{neuron_i}(:,4),comp_list{comp_i,2});

        if sum(valid_trials) < 5 || sum(nonvalid_trials) < 5
            nonviol_sdf_loop(comp_i,:) = nan(1,1001);
            viol_sdf_loop(comp_i,:) = nan(1,1001);
        else
            % Compute z-scored, smoothed SDFs for different identities at position 5
            viol_sdf_loop(comp_i,:) = [(smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials,1))),100)-baseline_mu_fr)./baseline_std_fr]';
            nonviol_sdf_loop(comp_i,:) = [(smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(nonvalid_trials,1))),100)-baseline_mu_fr)./baseline_std_fr]';
        end

        next_valid_trials = []; nonvalid_trials = [];
        % Identify valid trials: correct and labeled as 'nonviol'
        if comp_i == 1 || comp_i == 3
            next_valid_trials = find(cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'nonviol') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),comp_list{comp_i,1}) & strcmp(sdf_soundAlign_data{neuron_i}(:,4),comp_list{comp_i,2}))+1;
            next_nonvalid_trials = find(cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'viol') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),comp_list{comp_i,1}) & strcmp(sdf_soundAlign_data{neuron_i}(:,4),comp_list{comp_i,2}))+1;
        end

        if length(next_valid_trials) < 5 || length(next_nonvalid_trials) < 5 || any([isempty(next_valid_trials), isempty(next_nonvalid_trials)])
            next_viol_sdf_loop(comp_i,:) = nan(1,1001);
            next_nonviol_sdf_loop(comp_i,:) = nan(1,1001);
        else
            next_viol_sdf_loop(comp_i,:) = [(smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(next_valid_trials,1))),100)-baseline_mu_fr)./baseline_std_fr]';
            next_nonviol_sdf_loop(comp_i,:) = [(smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(next_nonvalid_trials,1))),100)-baseline_mu_fr)./baseline_std_fr]';
        end
    end

    viol_sdf(neuron_i,:) = nanmean(viol_sdf_loop);
    nonviol_sdf(neuron_i,:) = nanmean(nonviol_sdf_loop);

    next_viol_sdf(neuron_i,:) = nanmean(next_viol_sdf_loop);
    next_nonviol_sdf(neuron_i,:) = nanmean(next_nonviol_sdf_loop);

end

% Define time window of interest
time_win = [-100:563];

% Extract SDF signal for identity decoding in auditory and frontal areas
signal_in_violation_auditory = {viol_sdf(auditory_neuron_idx, 200+time_win),...
    nonviol_sdf(auditory_neuron_idx, 200+time_win)};
signal_in_violation_frontal = {viol_sdf(frontal_neuron_idx, 200+time_win),...
    nonviol_sdf(frontal_neuron_idx, 200+time_win)};

[pc_out_violation_auditory, ~] = get_xcond_pca(signal_in_violation_auditory);
[pc_out_violation_frontal, ~] = get_xcond_pca(signal_in_violation_frontal);


% Assign RGB colors for plotting
identity_color_1 = colorscale_identity(1,:); identity_color_2 = colorscale_identity(2,:); identity_color_3 = colorscale_identity(3,:);

% Plot 3D PCA trajectories
figuren('Renderer', 'painters', 'Position', [200,100,1200,400]);
subplot(1,2,1); hold on
plot3(pc_out_violation_auditory{1}(1,:)-pc_out_violation_auditory{1}(1,abs(time_win(1))),pc_out_violation_auditory{1}(2,:)-pc_out_violation_auditory{1}(2,abs(time_win(1))),pc_out_violation_auditory{1}(3,:)-pc_out_violation_auditory{1}(3,abs(time_win(1))),'LineWidth', plot_line_width, 'Color', identity_color_1)
plot3(pc_out_violation_auditory{2}(1,:)-pc_out_violation_auditory{2}(1,abs(time_win(1))),pc_out_violation_auditory{2}(2,:)-pc_out_violation_auditory{2}(2,abs(time_win(1))),pc_out_violation_auditory{2}(3,:)-pc_out_violation_auditory{2}(3,abs(time_win(1))),'LineWidth', plot_line_width, 'Color', identity_color_2)
scatter3(0,0,0,'ko','filled')
title('Violation (Aud)'); view(-45.1291,24.9773); grid on; legend({'violation','nonviolation'});

subplot(1,2,2); hold on
plot3(pc_out_violation_frontal{1}(1,:)-pc_out_violation_frontal{1}(1,abs(time_win(1))),pc_out_violation_frontal{1}(2,:)-pc_out_violation_frontal{1}(2,abs(time_win(1))),pc_out_violation_frontal{1}(3,:)-pc_out_violation_frontal{1}(3,abs(time_win(1))),'LineWidth', plot_line_width, 'Color', identity_color_1)
plot3(pc_out_violation_frontal{2}(1,:)-pc_out_violation_frontal{2}(1,abs(time_win(1))),pc_out_violation_frontal{2}(2,:)-pc_out_violation_frontal{2}(2,abs(time_win(1))),pc_out_violation_frontal{2}(3,:)-pc_out_violation_frontal{2}(3,abs(time_win(1))),'LineWidth', plot_line_width, 'Color', identity_color_2)
scatter3(0,0,0,'ko','filled')
title('Violation (Frontal)'); view(-45.1291,24.9773); grid on; legend({'violation','nonviolation'});


