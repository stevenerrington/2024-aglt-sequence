%% Does population activity at trial N's position-1 window decode trial N-1's group?
%
% prevtrial_confound_check.m established that trial N-1's sequence group
% predicts trial N's group (strong avoid-immediate-repeat structure:
% observed same-group repeat rate 0.161 vs 0.250 expected under
% independence, chi2(9)=77.2, p=6e-13) -- so a previous-trial carryover
% mechanism is statistically capable of explaining the identical-stimulus
% confound found in sequence_confound_check.m (frontal 89% / auditory 67%
% decoding of the CURRENT trial's group from position-1 activity, where
% the stimulus, sound "A", is physically identical across groups).
%
% This script tests the mechanism directly: does population activity at
% trial N's position-1 window (still the identical "A" stimulus) decode
% trial N-1's group? If it does, that is direct evidence of a genuine but
% mundane (non-anticipatory) carryover/adaptation effect from the
% previous trial, which combined with the previous-trial/current-trial
% statistical link above would be a complete, non-mysterious explanation
% for the original confound -- no prediction of the future required.
%
% Reuses the same extraction/decoding machinery as
% sequence_confound_check.m and prevtrial_confound_check.m. Assumes the
% same workspace variables: spike_log, dirs, auditory_neuron_idx,
% frontal_neuron_idx.

clear prevgroup_trials confound_results_prev

rng(20,'twister')

element_onset_ms = [0 563 1126 1688 2252];
time_win         = 0:413;
n_time_bins      = 4;
n_pcs            = 3;
min_trials       = 3;
n_pseudo_trial   = 15;
k_folds          = 5;
cv_repeats       = 4;
n_boot           = 300;
n_perm           = 300;

bin_edges = round(linspace(0, length(time_win), n_time_bins+1));

seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};
group_names     = {'prevgroup1','prevgroup2','prevgroup3','prevgroup4'};

n_neurons = size(spike_log,1);
prevgroup_trials(n_neurons,1) = struct();

for neuron_i = 1:n_neurons

    if mod(neuron_i,100) == 0
        fprintf('Neuron %i of %i\n', neuron_i, n_neurons);
    end

    try
        sdf_in = load(fullfile(dirs.root, 'data', 'spike', ...
            [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
        event_table_in = load(fullfile(dirs.mat_data, ...
            [spike_log.session{neuron_i} '.mat']), 'event_table');
    catch
        continue
    end

    event_table = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);

    trial_n_all   = event_table.trial_n;
    cond_value_all_full = event_table.cond_value;

    curr_group_full = nan(size(cond_value_all_full));
    for g = 1:4
        curr_group_full(ismember(cond_value_all_full, seq_cond_values{g})) = g;
    end
    curr_group_full(~nonviol_mask) = nan;

    % lookup: for literal trial_n value v, what group (if any) was it?
    max_trial_n = max(trial_n_all);
    group_by_trialn = nan(max_trial_n,1);
    group_by_trialn(trial_n_all) = curr_group_full;

    % previous-trial group for every row in this neuron's own trial table
    prev_group_full = nan(size(trial_n_all));
    has_prev = trial_n_all > 1;
    prev_group_full(has_prev) = group_by_trialn(max(trial_n_all(has_prev)-1,1));

    all_sdf = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    prev_group_nonviol = prev_group_full(nonviol_mask);

    baseline_window = 800:1000;
    baseline_mu_fr  = nanmean(nanmean(all_sdf(:,baseline_window)));
    baseline_std_fr = nanstd(nanmean(all_sdf(:,baseline_window)));

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue
    end

    for g = 1:4
        trial_idx = prev_group_nonviol == g;   % grouped by the PREVIOUS trial's identity
        if ~any(trial_idx)
            prevgroup_trials(neuron_i).(group_names{g}) = [];
            continue
        end
        cols = 1000 + element_onset_ms(1) + time_win;   % position 1 of the CURRENT trial -- identical "A" stimulus
        trial_window = all_sdf(trial_idx, cols);
        raw_bins = nan(size(trial_window,1), n_time_bins);
        for bin_i = 1:n_time_bins
            raw_bins(:,bin_i) = nanmean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2);
        end
        valid_rows = all(~isnan(raw_bins),2);
        raw_bins = raw_bins(valid_rows,:);
        if isempty(raw_bins)
            prevgroup_trials(neuron_i).(group_names{g}) = [];
        else
            prevgroup_trials(neuron_i).(group_names{g}) = (raw_bins - baseline_mu_fr) ./ baseline_std_fr;
        end
    end
end

%% Decode PREVIOUS trial's group (4-way) from CURRENT trial's identical position-1 window

decoding_problems = {
    'Prevgroup_from_position1_Auditory', auditory_neuron_idx;
    'Prevgroup_from_position1_Frontal',  frontal_neuron_idx;
    };

confound_results_prev = struct();

for prob_i = 1:size(decoding_problems,1)

    label      = decoding_problems{prob_i,1};
    neuron_idx = decoding_problems{prob_i,2};

    fprintf('\n=== %s ===\n', strrep(label,'_',' '));

    [boot_acc, perm_acc, n_used] = rigorous_lda_decode( ...
        prevgroup_trials, group_names, neuron_idx, min_trials, n_pseudo_trial, n_time_bins, ...
        n_pcs, k_folds, cv_repeats, n_boot, n_perm);

    confound_results_prev.(label).boot_acc     = boot_acc;
    confound_results_prev.(label).perm_acc     = perm_acc;
    confound_results_prev.(label).neurons_used = n_used;

    obs_acc = median(boot_acc,'omitnan');
    p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));

    fprintf('Neurons usable (>= %d real trials/class): %d / %d\n', min_trials, n_used, length(neuron_idx));
    fprintf('Nested-CV accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
    fprintf('Theoretical chance = 0.250 | Permutation null: median = %.3f, 95th pctile = %.3f\n', ...
        median(perm_acc,'omitnan'), prctile(perm_acc,95));
    fprintf('Permutation p-value (observed vs null): p = %.4f\n', p_perm);

end

fprintf('\nFor comparison, sequence_confound_check.m (decoding CURRENT trial''s group from the\n');
fprintf('same identical-stimulus window) found: Auditory 0.673, Frontal 0.892 (chance 0.250).\n');
fprintf('If THIS script''s numbers (decoding the PREVIOUS trial''s group) are comparably high,\n');
fprintf('that is direct evidence the original confound is carryover from trial N-1, not\n');
fprintf('anticipation of trial N.\n');


%% ================= Local functions (identical to the other confound-check scripts) =================

function [boot_acc, perm_acc, n_used] = rigorous_lda_decode(trial_struct, classes, neuron_idx, ...
    min_trials, n_pseudo_trial, n_time_bins, n_pcs, k_folds, cv_repeats, n_boot, n_perm)

n_classes = length(classes);

usable     = true(length(neuron_idx),1);
trial_pool = cell(length(neuron_idx), n_classes);

for ni = 1:length(neuron_idx)
    for ci = 1:n_classes
        v = trial_struct(neuron_idx(ni)).(classes{ci});
        trial_pool{ni,ci} = v;
        if isempty(v) || size(v,1) < min_trials
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
    y = y(randperm(length(y)));
    perm_acc(perm_i) = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats);
end

end


function [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins)

[n_neurons, n_classes] = size(trial_pool);
X = nan(n_pseudo_trial*n_classes, n_neurons*n_time_bins);
y = nan(n_pseudo_trial*n_classes, 1);

row0 = 0;
for ci = 1:n_classes
    block = nan(n_pseudo_trial, n_neurons*n_time_bins);
    for ni = 1:n_neurons
        pool = trial_pool{ni,ci};
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
