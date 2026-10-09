%% Confound check: can "which of the 4 sequences" be decoded from identical stimuli?
%
% Every one of the 4 grammatical sequences begins with the same physical
% sound, "A" (see sequence description.xlsx). So decoding "which sequence
% is this trial" from the position-1 window alone has NO possible acoustic
% driver -- there is nothing for the classifier to key off except
% non-specific, trial-level context (adaptation, drift, state). If this
% decodes well above chance, it demonstrates that the whole-sequence-
% anchored representation used in pca_lda_id_position_wholeseq_rigorous.m
% carries exactly this kind of confound, which would also inflate its
% "identity" decoding (C/F/G @ position 5 map almost 1:1 onto which
% sequence was played -- F only ever comes from sequence 2, G only from
% sequence 3, C only from sequences 1 or 4).
%
% Reuses the same extraction/decoding machinery as
% pca_lda_id_position_wholeseq_rigorous.m. Assumes the same workspace
% variables: spike_log, dirs, auditory_neuron_idx, frontal_neuron_idx.

clear seq_trials confound_results

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

seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};   % the 4 sequence groups (see wholeseq script header)
seq_group_names = {'seq1','seq2','seq3','seq4'};

n_neurons = size(spike_log,1);
seq_trials(n_neurons,1) = struct();

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

    baseline_window = 800:1000;
    all_sdf = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    cond_value_all = event_table.cond_value(nonviol_mask);

    baseline_mu_fr  = nanmean(nanmean(all_sdf(:,baseline_window)));
    baseline_std_fr = nanstd(nanmean(all_sdf(:,baseline_window)));

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue
    end

    for g = 1:4
        trial_idx = ismember(cond_value_all, seq_cond_values{g});
        if ~any(trial_idx)
            seq_trials(neuron_i).(seq_group_names{g}) = [];
            continue
        end
        cols = 1000 + element_onset_ms(1) + time_win;   % position 1 -- identical "A" stimulus every sequence
        trial_window = all_sdf(trial_idx, cols);
        raw_bins = nan(size(trial_window,1), n_time_bins);
        for bin_i = 1:n_time_bins
            raw_bins(:,bin_i) = nanmean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2);
        end
        valid_rows = all(~isnan(raw_bins),2);
        raw_bins = raw_bins(valid_rows,:);
        if isempty(raw_bins)
            seq_trials(neuron_i).(seq_group_names{g}) = [];
        else
            seq_trials(neuron_i).(seq_group_names{g}) = (raw_bins - baseline_mu_fr) ./ baseline_std_fr;
        end
    end
end

%% Decode sequence identity (4-way) from the identical position-1 "A" window, both areas

decoding_problems = {
    'Sequence_identity_Auditory', auditory_neuron_idx;
    'Sequence_identity_Frontal',  frontal_neuron_idx;
    };

confound_results = struct();

for prob_i = 1:size(decoding_problems,1)

    label      = decoding_problems{prob_i,1};
    neuron_idx = decoding_problems{prob_i,2};

    fprintf('\n=== %s (confound check: identical "A" stimulus) ===\n', strrep(label,'_',' '));

    [boot_acc, perm_acc, n_used] = rigorous_lda_decode( ...
        seq_trials, seq_group_names, neuron_idx, min_trials, n_pseudo_trial, n_time_bins, ...
        n_pcs, k_folds, cv_repeats, n_boot, n_perm);

    confound_results.(label).boot_acc     = boot_acc;
    confound_results.(label).perm_acc     = perm_acc;
    confound_results.(label).neurons_used = n_used;

    obs_acc = median(boot_acc,'omitnan');
    p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));

    fprintf('Neurons usable (>= %d real trials/class): %d / %d\n', min_trials, n_used, length(neuron_idx));
    fprintf('Nested-CV accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
    fprintf('Theoretical chance = 0.250 | Permutation null: median = %.3f, 95th pctile = %.3f\n', ...
        median(perm_acc,'omitnan'), prctile(perm_acc,95));
    fprintf('Permutation p-value (observed vs null): p = %.4f\n', p_perm);
    fprintf('--> If this is well above chance, population activity separates trials by sequence/context\n');
    fprintf('    even when the stimulus itself (sound "A") is physically identical.\n');

end


%% ================= Local functions (identical to pca_lda_id_position_wholeseq_rigorous.m) =================

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
