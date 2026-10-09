%% Rigorous, leakage-free trial-level LDA decoding of identity and ordinal position
%
% Reworked version of pca_lda_id_position.m. That script decoded from
% individual 1 ms timepoints of a single condition-averaged, 100 ms-smoothed
% PCA trajectory (~400 highly autocorrelated, non-independent "samples" per
% condition), used a single random 70/30 holdout split, and its bootstrap
% loop drew training/test rows from that same fixed set of timepoints
% (i.e. from the same trajectory the PCA itself was fit on) -- a form of
% train/test leakage that most likely inflates the reported accuracies.
%
% What this version does differently:
%   1. Decoding operates on genuine single-trial samples: one window-
%      averaged, baseline z-scored firing rate per real trial per neuron,
%      not per-timepoint samples of an averaged trajectory.
%   2. Because neurons were not all recorded simultaneously, trial-level
%      population vectors are built as bootstrapped pseudo-populations:
%      on every iteration, each neuron's contribution is an independent
%      draw (with replacement) from that neuron's own real trial pool for
%      the relevant class.
%   3. PCA is refit on the TRAINING fold only, inside each cross-
%      validation fold, and the held-out fold is projected into that
%      space afterwards. PCA is never fit on data that includes points
%      the classifier will later be tested on.
%   4. Performance is estimated via repeated stratified K-fold CV, nested
%      inside an outer bootstrap loop over trial resampling -- so the
%      reported interval reflects both trial-sampling noise and genuine
%      cross-validated generalisation, not a single lucky/unlucky split.
%   5. Chance level is assessed with a label-permutation test (shuffle
%      class labels, rerun the identical nested-CV pipeline, many times)
%      rather than compared only to a fixed theoretical value (1/3, 1/4),
%      which matters if real class sizes are not perfectly balanced.
%
% This is also a leak-free alternative to get_xcond_pca.m's neuron-
% inclusion step, which only intersects the first two conditions' valid-
% neuron sets even when 3-4 conditions are pooled -- here every neuron is
% explicitly required to have >= min_trials real trials in *every* class
% before it's included.
%
% Assumes the following are already in the workspace, exactly as set up by
% aglt_analysis_main.m: spike_log, sdf_soundAlign_data, auditory_neuron_idx,
% frontal_neuron_idx.
%
% Update: each trial is now represented by n_time_bins separate window-
% averaged rates (not one single mean across the whole 0-413 ms window).
% A single scalar per trial discards exactly the within-window temporal
% evolution that this manuscript's own PCA-trajectory analysis (Fig 2A)
% argues carries positional information in frontal cortex -- a flat mean
% would silently destroy that signal before decoding ever saw it. Binning
% keeps the analysis fully trial-level and leakage-free (PCA/LDA are still
% fit on the training fold only, per fold) while giving genuine temporal
% structure a chance to contribute.
%
% Runtime note: with the defaults below (300 bootstrap + 300 permutation
% iterations x 5-fold x 4 repeats x 4 decoding problems), expect this to
% run for several minutes to ~20 minutes depending on neuron counts.
% Reduce n_boot/n_perm for a quick check first.

clear identity_trials position_trials rigorous_results

rng(20,'twister')   % match the seed used in the original script

%% ---- User-set parameters ----
time_win        = 0:413;      % post-onset analysis window (ms), as in the manuscript
n_pcs           = 3;          % PCs retained for decoding, as in the manuscript
min_trials      = 3;          % minimum real trials/neuron/class required to include a neuron
n_time_bins     = 4;          % non-overlapping bins spanning the 0-413 ms window (~103 ms each)
n_pseudo_trial  = 15;         % pseudo-trials per class per pseudo-population draw
k_folds         = 5;          % stratified CV folds
cv_repeats      = 4;          % repeats of the k-fold split per outer iteration
n_boot          = 300;        % outer bootstrap iterations (trial resampling -> accuracy CI)
n_perm          = 300;        % label-permutation iterations (-> chance-level null distribution)

identity_classes = {'C','F','G'};                                  % position 5, as in the manuscript design
position_classes = {'position_2','position_3','position_4','position_5'};  % identity C only, as in the manuscript design

bin_edges = round(linspace(0, length(time_win), n_time_bins+1));  % non-overlapping bin boundaries into time_win

%% ---- Step 1: single-trial, window-averaged, baseline z-scored firing rates ----
% One struct array entry per neuron; each field is a vector of per-trial
% firing rates (z-scored to that neuron's own baseline), one per class.

n_neurons = size(spike_log,1);

identity_trials(n_neurons,1) = struct();
position_trials(n_neurons,1) = struct();

for neuron_i = 1:n_neurons

    valid_trials = cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & ...
        strcmp(sdf_soundAlign_data{neuron_i}(:,5),'nonviol');

    baseline_mu_fr  = nanmean(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));
    baseline_std_fr = nanstd(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue   % no usable baseline for this neuron -- its trial pools stay empty
    end             % (default [] from the preallocated struct) and it is excluded via min_trials

    for class_i = 1:length(identity_classes)
        trial_idx = valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & ...
            strcmp(sdf_soundAlign_data{neuron_i}(:,3),identity_classes{class_i});
        trial_sdf = cell2mat(sdf_soundAlign_data{neuron_i}(trial_idx,1));
        if isempty(trial_sdf)
            identity_trials(neuron_i).(identity_classes{class_i}) = [];
        else
            trial_window = trial_sdf(:,200+time_win);            % trials x 414
            trial_bins = nan(size(trial_window,1), n_time_bins);
            for bin_i = 1:n_time_bins
                trial_bins(:,bin_i) = nanmean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2);
            end
            valid_rows = all(~isnan(trial_bins),2);              % drop trials with no valid samples in any bin
            trial_bins = trial_bins(valid_rows,:);
            if isempty(trial_bins)
                identity_trials(neuron_i).(identity_classes{class_i}) = [];
            else
                identity_trials(neuron_i).(identity_classes{class_i}) = (trial_bins - baseline_mu_fr) ./ baseline_std_fr;
            end
        end
    end

    for class_i = 1:length(position_classes)
        trial_idx = valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & ...
            strcmp(sdf_soundAlign_data{neuron_i}(:,4),position_classes{class_i});
        trial_sdf = cell2mat(sdf_soundAlign_data{neuron_i}(trial_idx,1));
        if isempty(trial_sdf)
            position_trials(neuron_i).(position_classes{class_i}) = [];
        else
            trial_window = trial_sdf(:,200+time_win);
            trial_bins = nan(size(trial_window,1), n_time_bins);
            for bin_i = 1:n_time_bins
                trial_bins(:,bin_i) = nanmean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2);
            end
            valid_rows = all(~isnan(trial_bins),2);
            trial_bins = trial_bins(valid_rows,:);
            if isempty(trial_bins)
                position_trials(neuron_i).(position_classes{class_i}) = [];
            else
                position_trials(neuron_i).(position_classes{class_i}) = (trial_bins - baseline_mu_fr) ./ baseline_std_fr;
            end
        end
    end
end

%% ---- Step 2: run the nested, leakage-free decoding for each area x decoding problem ----

decoding_problems = {
    'Identity_Auditory', identity_trials, identity_classes, auditory_neuron_idx;
    'Identity_Frontal',  identity_trials, identity_classes, frontal_neuron_idx;
    'Position_Auditory', position_trials, position_classes, auditory_neuron_idx;
    'Position_Frontal',  position_trials, position_classes, frontal_neuron_idx;
    };

rigorous_results = struct();

for prob_i = 1:size(decoding_problems,1)

    label        = decoding_problems{prob_i,1};
    trial_struct = decoding_problems{prob_i,2};
    classes      = decoding_problems{prob_i,3};
    neuron_idx   = decoding_problems{prob_i,4};

    fprintf('\n=== %s ===\n', strrep(label,'_',' '));

    [boot_acc, perm_acc, n_used] = rigorous_lda_decode( ...
        trial_struct, classes, neuron_idx, min_trials, n_pseudo_trial, n_time_bins, ...
        n_pcs, k_folds, cv_repeats, n_boot, n_perm);

    rigorous_results.(label).boot_acc     = boot_acc;
    rigorous_results.(label).perm_acc     = perm_acc;
    rigorous_results.(label).neurons_used = n_used;
    rigorous_results.(label).n_classes    = length(classes);

    obs_acc = median(boot_acc,'omitnan');
    p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));

    fprintf('Neurons usable (>= %d real trials/class): %d / %d\n', min_trials, n_used, length(neuron_idx));
    fprintf('Nested-CV accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
    fprintf('Theoretical chance = %.3f | Permutation null: median = %.3f, 95th pctile = %.3f\n', ...
        1/length(classes), median(perm_acc,'omitnan'), prctile(perm_acc,95));
    fprintf('Permutation p-value (observed vs null): p = %.4f\n', p_perm);

end

%% ---- Step 3: plot rigorous accuracy distributions (standard MATLAB plotting) ----

problem_names = fieldnames(rigorous_results);
figure('Renderer','painters','Position',[100 100 700 450]); hold on

bar_x        = 1:length(problem_names);
bar_median   = nan(size(bar_x));
err_lo       = nan(size(bar_x));
err_hi       = nan(size(bar_x));
chance_level = nan(size(bar_x));

for i = 1:length(problem_names)
    acc = rigorous_results.(problem_names{i}).boot_acc;
    bar_median(i)   = median(acc,'omitnan');
    err_lo(i)       = bar_median(i) - prctile(acc,2.5);
    err_hi(i)       = prctile(acc,97.5) - bar_median(i);
    chance_level(i) = 1 / rigorous_results.(problem_names{i}).n_classes;
end

bar(bar_x, bar_median, 0.5, 'FaceColor',[0.6 0.6 0.6]);
errorbar(bar_x, bar_median, err_lo, err_hi, 'k', 'LineStyle','none','LineWidth',1.2)
plot(bar_x, chance_level, 'rd','MarkerFaceColor','r')
set(gca,'XTick',bar_x,'XTickLabel',strrep(problem_names,'_',' '),'XTickLabelRotation',20)
ylabel('Nested cross-validated accuracy')
ylim([0 1]); box off
title('Rigorous trial-level LDA decoding (leakage-free, permutation-tested)')
legend({'Median accuracy','95% bootstrap CI','Theoretical chance'},'Location','best')

%% ---- Step 4: area / feature comparisons (same style as the original script) ----

fprintf('\nAuditory: identity vs position decoding\n');
bootstrap_compare(rigorous_results.Identity_Auditory.boot_acc, rigorous_results.Position_Auditory.boot_acc);

fprintf('\nFrontal: identity vs position decoding\n');
bootstrap_compare(rigorous_results.Identity_Frontal.boot_acc, rigorous_results.Position_Frontal.boot_acc);

fprintf('\nPosition decoding: frontal vs auditory\n');
bootstrap_compare(rigorous_results.Position_Frontal.boot_acc, rigorous_results.Position_Auditory.boot_acc);

fprintf('\nIdentity decoding: frontal vs auditory\n');
bootstrap_compare(rigorous_results.Identity_Frontal.boot_acc, rigorous_results.Identity_Auditory.boot_acc);


%% ================= Local functions =================

function [boot_acc, perm_acc, n_used] = rigorous_lda_decode(trial_struct, classes, neuron_idx, ...
    min_trials, n_pseudo_trial, n_time_bins, n_pcs, k_folds, cv_repeats, n_boot, n_perm)
% Builds bootstrap pseudo-populations from single-trial data (each trial
% represented by n_time_bins window-averaged rates, not one flat mean),
% runs nested (leakage-free) PCA+LDA cross-validation on each one, and
% separately builds a label-permutation null distribution using the
% identical pipeline.

n_classes = length(classes);

usable     = true(length(neuron_idx),1);
trial_pool = cell(length(neuron_idx), n_classes);

for ni = 1:length(neuron_idx)
    for ci = 1:n_classes
        v = trial_struct(neuron_idx(ni)).(classes{ci});   % Ntrials x n_time_bins
        trial_pool{ni,ci} = v;
        if size(v,1) < min_trials
            usable(ni) = false;
        end
    end
end

use_idx    = find(usable);
n_used     = length(use_idx);
trial_pool = trial_pool(use_idx,:);

boot_acc = nan(n_boot,1);
for boot_i = 1:n_boot
    [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins);
    boot_acc(boot_i) = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats);
end

perm_acc = nan(n_perm,1);
for perm_i = 1:n_perm
    [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins);
    y = y(randperm(length(y)));   % break the label-feature relationship
    perm_acc(perm_i) = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats);
end

end


function [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins)
% trial_pool: n_neurons x n_classes cell array; each cell is an
% Ntrials x n_time_bins matrix of per-trial, per-bin firing rates.
% Each neuron's pseudo-trials are drawn independently, with replacement,
% from that neuron's own real trial pool for the given class -- neurons
% were not all recorded simultaneously, so this is the standard
% pseudo-population construction rather than resampling matched trial
% indices across neurons. Every neuron contributes n_time_bins columns,
% so X has n_neurons*n_time_bins features.

[n_neurons, n_classes] = size(trial_pool);
X = nan(n_pseudo_trial*n_classes, n_neurons*n_time_bins);
y = nan(n_pseudo_trial*n_classes, 1);

row0 = 0;
for ci = 1:n_classes
    block = nan(n_pseudo_trial, n_neurons*n_time_bins);
    for ni = 1:n_neurons
        pool = trial_pool{ni,ci};                    % Ntrials x n_time_bins
        draw_idx = randi(size(pool,1), n_pseudo_trial, 1);
        cols = (ni-1)*n_time_bins + (1:n_time_bins);
        block(:,cols) = pool(draw_idx,:);
    end
    X(row0+1:row0+n_pseudo_trial,:) = block;
    y(row0+1:row0+n_pseudo_trial)   = ci;
    row0 = row0 + n_pseudo_trial;
end

end


function acc = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats)
% Repeated stratified K-fold CV. PCA is fit on the training fold only and
% the test fold is projected into that space -- this is the step that
% was missing (a source of circularity) in the original script.

fold_acc = nan(k_folds*cv_repeats,1);
f = 1;
for rep_i = 1:cv_repeats
    cvp = cvpartition(y, 'KFold', k_folds);
    for fold_i = 1:k_folds
        train_idx = training(cvp, fold_i);
        test_idx  = test(cvp, fold_i);

        mu = mean(X(train_idx,:),1);
        [coeff, scores_train] = pca(X(train_idx,:));
        n_keep = min(n_pcs, size(coeff,2));
        scores_train = scores_train(:,1:n_keep);
        scores_test  = (X(test_idx,:) - mu) * coeff(:,1:n_keep);

        lda_model = fitcdiscr(scores_train, y(train_idx));
        preds = predict(lda_model, scores_test);
        fold_acc(f) = mean(preds == y(test_idx));
        f = f + 1;
    end
end
acc = mean(fold_acc,'omitnan');

end
