% Clear previously stored SDF-related variables to avoid conflicts
clear a_sdf c_sdf d_sdf f_sdf g_sdf pos1_sdf pos2_sdf pos3_sdf pos4_sdf pos5_sdf

% Loop over neurons
for neuron_i = 1:size(spike_log,1)

    % Identify valid trials:
    % Column 9 = correct trials (logical), Column 5 = 'nonviol'
    valid_trials = cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'nonviol');

    % Compute baseline firing activity (mean and SD) for z-scoring
    baseline_mu_fr = nanmean(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));
    baseline_std_fr = nanstd(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));

    % ============================
    % Compute z-scored SDFs for identity decoding at POSITION 5
    % (identities A, C, D, F, G)
    % ============================
    a_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    c_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    d_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    f_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    g_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    % ============================
    % Compute z-scored SDFs for POSITION decoding for identity C
    % (positions 1–5)
    % ============================
    pos1_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    pos2_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    pos3_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    pos4_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    pos5_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    % ============================
    % Count number of trials for identity decoding at POSITION 4
    % ============================
    n_trl_id(neuron_i,1) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'));
    n_trl_id(neuron_i,2) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_id(neuron_i,3) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'));
    n_trl_id(neuron_i,4) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'));
    n_trl_id(neuron_i,5) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'));

    % ============================
    % Count number of trials for POSITION decoding for identity C
    % ============================
    n_trl_pos(neuron_i,1) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,2) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,3) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,4) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,5) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));

end

% Time window used for PCA extraction (relative to alignment)
time_win = 0:413;

% Extract identity SDFs for auditory neurons (C, F, G at position 5)
signal_in_identity_auditory = {c_sdf(auditory_neuron_idx, 200+time_win),...
    f_sdf(auditory_neuron_idx, 200+time_win),...
    g_sdf(auditory_neuron_idx, 200+time_win)};

% Same for frontal neurons
signal_in_identity_frontal = {c_sdf(frontal_neuron_idx, 200+time_win),...
    f_sdf(frontal_neuron_idx, 200+time_win),...
    g_sdf(frontal_neuron_idx, 200+time_win)};

% Extract position decoding SDFs (positions 2–5)
signal_in_order_auditory = {pos2_sdf(auditory_neuron_idx, 200+time_win),...
    pos3_sdf(auditory_neuron_idx, 200+time_win),...
    pos4_sdf(auditory_neuron_idx, 200+time_win),...
    pos5_sdf(auditory_neuron_idx, 200+time_win)};

signal_in_order_frontal = {pos2_sdf(frontal_neuron_idx, 200+time_win),...
    pos3_sdf(frontal_neuron_idx, 200+time_win),...
    pos4_sdf(frontal_neuron_idx, 200+time_win),...
    pos5_sdf(frontal_neuron_idx, 200+time_win)};

% Perform PCA for each group of conditions
[pc_out_identity_auditory, ~] = get_xcond_pca(signal_in_identity_auditory);
[pc_out_identity_frontal, ~] = get_xcond_pca(signal_in_identity_frontal);
[pc_out_order_auditory, ~] = get_xcond_pca(signal_in_order_auditory);
[pc_out_order_frontal, ~] = get_xcond_pca(signal_in_order_frontal);

% Fix random seed for reproducible cross-validation and bootstrapping
rng(20)

% Store PCA outputs for looping
pc_list_in = {pc_out_identity_auditory, pc_out_identity_frontal, pc_out_order_auditory, pc_out_order_frontal};

pc_to_include = {[1], [1,2,3], [2], [2 3], [2 3 4], [2 3 4 5], [2 3 4 5 6], [2 3 4 5 6 7], [2 3 4 5 6 7 8], [2 3 4 5 6 7 8 9], [2 3 4 5 6 7 8 9 10]};

% ============================
% LDA CLASSIFICATION
% ============================
for data_i = 1:length(pc_list_in)

    fprintf('Extracting data from data_i = %i \n', data_i)
    % Load PCA result for this dataset
    pc_in = pc_list_in{data_i};
    nconds = size(pc_in,2);
    npcs = 3;  % use first 3 PCs

    for pc_cond_i = 1:length(pc_to_include)
        fprintf('       Extracting data from pc_cond_i = %i \n', pc_cond_i)

        % Flatten PCA trajectories → rows = timepoints, columns = PCs
        scores = []; labels = [];
        for cond_i = 1:nconds
            scores_a = [];
            for pc_i = pc_to_include{pc_cond_i}
                scores_a = [scores_a, pc_in{cond_i}(pc_i,:)'];  % combine PCs
            end
            scores = [scores; scores_a];                      % append
            labels = [labels, repmat(cond_i,1,length(scores_a))]; % condition labels
        end

        % Hold-out CV: train 70%, test 30%
        cv = cvpartition(labels, 'Holdout', 0.3);
        train_idx = training(cv);
        test_idx = test(cv);

        % Fit LDA model
        lda_model = fitcdiscr(scores(train_idx,:), labels(train_idx));

        % Predict and compute accuracy
        preds = predict(lda_model, scores(test_idx,:));
        cv_accuracy = sum(preds' == labels(test_idx)) / sum(test_idx);
        cv_accuracy_out(data_i) = cv_accuracy;

        % Normalized confusion matrix
        cm = confusionmat(labels(test_idx), preds);
        cm_normalized = cm ./ sum(cm, 2);
        classifier_confusionMatrix{data_i} = cm_normalized;

        % ============================
        % Bootstrapped classifier accuracy
        % ============================
        nBoots = 1000;
        acc_boot = nan(nBoots, 1);

        for boot_i = 1:nBoots
            n = size(scores,1);

            % Resample with replacement for training
            train_idx = randsample(n, round(0.7*n), true);
            test_idx = setdiff(1:n, unique(train_idx));

            % Skip degenerate samples
            if isempty(test_idx) || length(unique(labels(train_idx))) < 2
                continue
            end

            lda_model = fitcdiscr(scores(train_idx,:), labels(train_idx));
            predicted = predict(lda_model, scores(test_idx,:));
            acc_boot(boot_i) = mean(predicted == labels(test_idx)');
        end

        classifier_bootstrapAccuracy{data_i, pc_cond_i} = acc_boot;
    end
end


%%
figure('Renderer', 'painters', 'Position', [200 200 1200 250]); hold on;
nData = length(pc_list_in);
nCond = length(pc_to_include);

pc_list_in = {pc_out_identity_auditory, pc_out_identity_frontal, pc_out_order_auditory, pc_out_order_frontal};

cond_labels = {'Aud - ID', 'Frontal - ID', 'Auditory - Order','Frontal - Order'};

subplot_count = 0;

for data_i = [4 2 3 1]
    subplot_count = subplot_count + 1;
    subplot(1, nData, subplot_count); hold on;

    means = zeros(1, nCond);
    ci_low = zeros(1, nCond);
    ci_high = zeros(1, nCond);

    for pc_cond_i = 1:nCond
        vals = classifier_bootstrapAccuracy{data_i, pc_cond_i};

        means(pc_cond_i) = mean(vals);

        ci = prctile(vals, [2.5 97.5]); % 95% CI
        ci_low(pc_cond_i) = ci(1);
        ci_high(pc_cond_i) = ci(2);
    end

    x = 1:nCond;

    % Mean line
    plot(x, means, '-o', 'LineWidth', 1.5);

    % Shaded CI
    fill([x fliplr(x)], ...
         [ci_low fliplr(ci_high)], ...
         [0.8 0.8 1], ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');

    xlabel('Number of included PCs');
    ylabel('Accuracy');
    title(cond_labels{data_i});
    xticks(x);
    ylim([0.1 1.0])
    grid on;
end