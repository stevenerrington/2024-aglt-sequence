% Clear previously stored data variables related to SDFs
clear a_sdf c_sdf d_sdf f_sdf g_sdf pos1_sdf pos2_sdf pos3_sdf pos4_sdf pos5_sdf

% Loop over each neuron in spike_log
for neuron_i = 1:size(spike_log,1)

    % Identify valid trials: correct and labeled as 'nonviol'
    valid_trials = cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'nonviol');

    % Compute baseline firing rate (mean and std) for z-scoring
    baseline_mu_fr = nanmean(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));
    baseline_std_fr = nanstd(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));

    % Compute z-scored, smoothed SDFs for different identities at position 5
    a_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    c_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    d_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    f_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    g_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    % Compute z-scored, smoothed SDFs for identity 'C' across all positions
    pos1_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos2_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos3_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos4_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos5_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    % Count valid trials for each identity at position 4
    n_trl_id(neuron_i,1) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'));
    n_trl_id(neuron_i,2) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_id(neuron_i,3) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'));
    n_trl_id(neuron_i,4) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'));
    n_trl_id(neuron_i,5) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'));

    % Count valid trials for identity 'C' across all positions
    n_trl_pos(neuron_i,1) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,2) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,3) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,4) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,5) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));

end

% Define time window of interest
time_win = [0:413];

% Extract SDF signal for identity decoding in auditory and frontal areas
signal_in_identity_auditory = {c_sdf(auditory_neuron_idx, 200+time_win),...
    f_sdf(auditory_neuron_idx, 200+time_win),...
    g_sdf(auditory_neuron_idx, 200+time_win)};

signal_in_identity_frontal = {c_sdf(frontal_neuron_idx, 200+time_win),...
    f_sdf(frontal_neuron_idx, 200+time_win),...
    g_sdf(frontal_neuron_idx, 200+time_win)};

% Extract SDF signal for order decoding in auditory and frontal areas
signal_in_order_auditory = {pos2_sdf(auditory_neuron_idx, 200+time_win),...
    pos3_sdf(auditory_neuron_idx, 200+time_win),...
    pos4_sdf(auditory_neuron_idx, 200+time_win),...
    pos5_sdf(auditory_neuron_idx, 200+time_win)};

signal_in_order_frontal = {pos2_sdf(frontal_neuron_idx, 200+time_win),...
    pos3_sdf(frontal_neuron_idx, 200+time_win),...
    pos4_sdf(frontal_neuron_idx, 200+time_win),...
    pos5_sdf(frontal_neuron_idx, 200+time_win)};

% Run PCA on all sets of data
[pc_out_identity_auditory, ~] = get_xcond_pca(signal_in_identity_auditory);
[pc_out_identity_frontal, ~] = get_xcond_pca(signal_in_identity_frontal);
[pc_out_order_auditory, ~] = get_xcond_pca(signal_in_order_auditory);
[pc_out_order_frontal, ~] = get_xcond_pca(signal_in_order_frontal);

% Initialize RNG for reproducibility
rng(20)

% Prepare list of PCA outputs
clear pc_list_in
pc_list_in = {pc_out_identity_auditory, pc_out_identity_frontal, pc_out_order_auditory, pc_out_order_frontal};

% Loop through each PCA dataset for classification
for data_i = 1:length(pc_list_in)

    % Extract PCs for current dataset
    clear pc_in
    pc_in = pc_list_in{data_i};
    nconds = size(pc_in,2); npcs = 3;

    % Restructure PCA scores into feature matrix and labels
    scores = []; labels = [];
    for cond_i = 1:nconds
        scores_a = [];
        for pc_i = 1:npcs
            scores_a = [scores_a, pc_in{cond_i}(pc_i,:)'];
        end
        scores = [scores; scores_a];
        labels = [labels, repmat(cond_i,1,length(scores_a))];
    end

    % Set up and train LDA classifier with cross-validation
    clear cv cv_accuracy all_preds lda_model preds
    cv = cvpartition(labels, 'Holdout', 0.3);
    train_idx = training(cv); test_idx = test(cv);

    lda_model = fitcdiscr(scores(train_idx,:), labels(train_idx));
    preds = predict(lda_model, scores(test_idx,:));
    cv_accuracy = sum(preds' == labels(test_idx)) / sum(test_idx);

    % Normalize confusion matrix
    cm = confusionmat(labels(test_idx), preds);
    cm_normalized = cm ./ sum(cm, 2);
    classifier_confusionMatrix{data_i} = cm_normalized;

    % Run bootstrap to estimate classifier performance
    nBoots = 1000;
    acc_boot = nan(nBoots, 1);

    for boot_i = 1:nBoots
        n = size(scores,1);
        train_idx = randsample(n, round(0.7*n), true);
        test_idx = setdiff(1:n, unique(train_idx));

        if isempty(test_idx) || length(unique(labels(train_idx))) < 2
            continue
        end

        lda_model = fitcdiscr(scores(train_idx,:), labels(train_idx));
        predicted = predict(lda_model, scores(test_idx,:));
        acc_boot(boot_i) = mean(predicted == labels(test_idx)');
    end

    classifier_bootstrapAccuracy{data_i} = acc_boot;
end

% Visualization setup for PCA and confusion matrices
plot_line_width = 1;
colorscale_identity = cbrewer('qual','Set1',3);
colorscale_position = cbrewer('qual','Set2',4);
colorscale_heatmap = cbrewer('seq','Greens',100);
colorscale_heatmap(colorscale_heatmap < 0) = 0;

% Assign RGB colors for plotting
identity_color_1 = colorscale_identity(1,:); identity_color_2 = colorscale_identity(2,:); identity_color_3 = colorscale_identity(3,:);
position_color_1 = colorscale_position(1,:); position_color_2 = colorscale_position(2,:); position_color_3 = colorscale_position(3,:); position_color_4 = colorscale_position(4,:);

% Plot 3D PCA trajectories
figuren('Renderer', 'painters', 'Position', [200,100,1200,683]);
subplot(2,4,1); hold on
plot3(pc_out_identity_auditory{1}(1,:),pc_out_identity_auditory{1}(2,:),pc_out_identity_auditory{1}(3,:),'LineWidth', plot_line_width, 'Color', identity_color_1)
plot3(pc_out_identity_auditory{2}(1,:),pc_out_identity_auditory{2}(2,:),pc_out_identity_auditory{2}(3,:),'LineWidth', plot_line_width, 'Color', identity_color_2)
plot3(pc_out_identity_auditory{3}(1,:),pc_out_identity_auditory{3}(2,:),pc_out_identity_auditory{3}(3,:),'LineWidth', plot_line_width, 'Color', identity_color_3)
title('Identity (Aud)'); view(-45.1291,24.9773); grid on

subplot(2,4,2); hold on
plot3(pc_out_identity_frontal{1}(1,:),pc_out_identity_frontal{1}(2,:),pc_out_identity_frontal{1}(3,:),'LineWidth', plot_line_width, 'Color', identity_color_1)
plot3(pc_out_identity_frontal{2}(1,:),pc_out_identity_frontal{2}(2,:),pc_out_identity_frontal{2}(3,:),'LineWidth', plot_line_width, 'Color', identity_color_2)
plot3(pc_out_identity_frontal{3}(1,:),pc_out_identity_frontal{3}(2,:),pc_out_identity_frontal{3}(3,:),'LineWidth', plot_line_width, 'Color', identity_color_3)
title('Identity (Frontal)'); view(-45.1291,24.9773); grid on

subplot(2,4,3); hold on
plot3(pc_out_order_auditory{1}(1,:),pc_out_order_auditory{1}(2,:),pc_out_order_auditory{1}(3,:),'LineWidth', plot_line_width, 'Color', position_color_1)
plot3(pc_out_order_auditory{2}(1,:),pc_out_order_auditory{2}(2,:),pc_out_order_auditory{2}(3,:),'LineWidth', plot_line_width, 'Color', position_color_2)
plot3(pc_out_order_auditory{3}(1,:),pc_out_order_auditory{3}(2,:),pc_out_order_auditory{3}(3,:),'LineWidth', plot_line_width, 'Color', position_color_3)
plot3(pc_out_order_auditory{4}(1,:),pc_out_order_auditory{4}(2,:),pc_out_order_auditory{4}(3,:),'LineWidth', plot_line_width, 'Color',position_color_4)
title('Position (Aud)'); view(-45.1291,24.9773); grid on

subplot(2,4,4); hold on
plot3(pc_out_order_frontal{1}(1,:),pc_out_order_frontal{1}(2,:),pc_out_order_frontal{1}(3,:),'LineWidth', plot_line_width, 'Color', position_color_1)
plot3(pc_out_order_frontal{2}(1,:),pc_out_order_frontal{2}(2,:),pc_out_order_frontal{2}(3,:),'LineWidth', plot_line_width, 'Color', position_color_2)
plot3(pc_out_order_frontal{3}(1,:),pc_out_order_frontal{3}(2,:),pc_out_order_frontal{3}(3,:),'LineWidth', plot_line_width, 'Color', position_color_3)
plot3(pc_out_order_frontal{4}(1,:),pc_out_order_frontal{4}(2,:),pc_out_order_frontal{4}(3,:),'LineWidth', plot_line_width, 'Color', position_color_4)
title('Position (Frontal)'); view(-45.1291,24.9773); grid on


% Plot confusion matrices
for data_i = 1:4
    subplot(2,4,data_i+4)
    imagesc(classifier_confusionMatrix{data_i})
    grid off; box off
    colorbar('SouthOutside')  % Places the colorbar below the matrix
    clim([0 1]); colormap(colorscale_heatmap)
    axis square
end

% Aggregate bootstrap results for bar plot
plot_bootstrap_data = []; plot_bootstrap_cond = []; plot_bootstrap_area = [];

for data_i = 1:4
    plot_bootstrap_data = [plot_bootstrap_data; classifier_bootstrapAccuracy{data_i}];

    switch data_i
        case 1; cond = 'identity'; area = 'auditory';
        case 2; cond = 'identity'; area = 'frontal';
        case 3; cond = 'order';    area = 'auditory';
        case 4; cond = 'order';    area = 'frontal';
    end

    plot_bootstrap_cond = [plot_bootstrap_cond; repmat({cond},length(classifier_bootstrapAccuracy{data_i}),1)];
    plot_bootstrap_area = [plot_bootstrap_area; repmat({area},length(classifier_bootstrapAccuracy{data_i}),1)];
end

% Plot bootstrap classification accuracy
figure('Renderer', 'painters', 'Position', [100 100 500 500]);
plot_bootstrap_acc(1,1) = gramm('x', plot_bootstrap_area, 'y', plot_bootstrap_data, 'color', plot_bootstrap_cond);
plot_bootstrap_acc(1,1).stat_summary('geom', {'bar', 'black_errorbar'}, 'dodge', 1, 'width', 0.5);
plot_bootstrap_acc(1,1).geom_hline('yintercept',0.25,'style','-.');
plot_bootstrap_acc(1,1).geom_hline('yintercept',0.33,'style','--');
plot_bootstrap_acc(1,1).axe_property('YLim', [0 1]);
plot_bootstrap_acc.draw;
