%% Cross-sequence-group generalization test of ordinal position decoding
%
% prevtrial_neural_check.m established that population activity carries a
% strong, long-lived (trial-N-1-carryover) sequence-CONTEXT signal, much
% stronger in frontal cortex than auditory cortex, that is entirely
% independent of the physical stimulus. That is a serious problem for the
% manuscript's position-decoding claim specifically -- not just the
% confound-check -- because with only 4 fixed grammatical sequences,
% "which ordinal position is letter C at" and "which sequence is this"
% are NOT independent:
%
%   seq1 (cond_value 1,5) A C G F C   -- C at position 2, 5
%   seq2 (cond_value 2,6) A D C G F   -- C at position 3
%   seq3 (cond_value 3,7) A C F C G   -- C at position 2, 4
%   seq4 (cond_value 4,8) A D C F C   -- C at position 3, 5
%
% Position 4-C occurs ONLY in seq3 -- a full, structural confound with
% sequence group that no analysis can undo, so it is dropped here.
% Positions 2, 3 and 5 each span exactly two sequences, which makes a
% clean, structural control possible: split the 4 sequences into two
% group-sets that each cover all three positions using DIFFERENT
% physical sequences --
%   Set A = sequences {1, 2}   (supplies position 2 via seq1, position 3
%                                via seq2, position 5 via seq1)
%   Set B = sequences {3, 4}   (supplies position 2 via seq3, position 3
%                                via seq4, position 5 via seq4)
% -- then train the classifier on Set A only and test on Set B only (and
% vice versa). Because train and test trials come from entirely
% different physical sequences, any sequence-group-identity signal
% (including the carryover contamination) cannot transfer across the
% split and cannot be used as a decoding shortcut. Only genuine,
% sequence-independent ordinal-position coding can generalize.
%
% (Identity decoding, C vs F vs G at position 5, cannot be tested this
% way: F only ever occurs in seq2 and G only ever occurs in seq3, so 2 of
% 3 classes have no cross-group alternative at all -- that analysis is
% out of scope here.)
%
% Reuses the same continuous, sequence-onset-anchored feature extraction
% as pca_lda_id_position_wholeseq_rigorous.m. Assumes the following are
% already in the workspace: spike_log, dirs, auditory_neuron_idx,
% frontal_neuron_idx.

clear trials_A trials_B crossgroup_results

rng(20,'twister')

%% ---- User-set parameters ----
element_onset_ms = [0 563 1126 1688 2252];
time_win         = 0:413;
n_time_bins      = 4;
n_pcs            = 3;
min_trials       = 3;
n_pseudo_trial   = 15;
n_boot           = 300;
n_perm           = 300;

bin_edges = round(linspace(0, length(time_win), n_time_bins+1));

seq_letters = { ...
    'A','C','G','F','C'; ...   % cond_value 1 & 5 -- Grammatical5_novel   (Set A)
    'A','D','C','G','F'; ...   % cond_value 2 & 6 -- Grammatical10_novel (Set A)
    'A','C','F','C','G'; ...   % cond_value 3 & 7 -- Grammatical_3       (Set B)
    'A','D','C','F','C'};      % cond_value 4 & 8 -- Grammatical_8       (Set B)
seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};

groupset_A = [1 2];
groupset_B = [3 4];

position_list     = [2 3 5];   % position 4 dropped -- structurally confounded with seq3, see header
position_classes  = {'position_2','position_3','position_5'};

% Sanity-print which physical sequence feeds each class in each set, so
% the group assignment above can be checked at a glance.
fprintf('Cross-group assignment (position -> group used in Set A / Set B):\n');
for i = 1:length(position_list)
    pos = position_list(i);
    cand = find(strcmp(seq_letters(:,pos), 'C'))';
    gA = intersect(cand, groupset_A);
    gB = intersect(cand, groupset_B);
    fprintf('  position %d : Set A = seq%d, Set B = seq%d\n', pos, gA(1), gB(1));
end

%% ---- Step 1: single-trial, whole-sequence-anchored, baseline z-scored bin features ----
% Two parallel struct arrays: trials_A drawn only from Set-A sequences,
% trials_B only from Set-B sequences -- disjoint physical trials, never
% mixed within a class.

n_neurons = size(spike_log,1);
trials_A(n_neurons,1) = struct();
trials_B(n_neurons,1) = struct();

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

    for i = 1:length(position_list)
        pos = position_list(i);
        cond_values_A = lookup_cond_values_groupset(seq_letters, seq_cond_values, pos, 'C', groupset_A);
        cond_values_B = lookup_cond_values_groupset(seq_letters, seq_cond_values, pos, 'C', groupset_B);

        trials_A(neuron_i).(position_classes{i}) = extract_bin_features( ...
            all_sdf, cond_value_all, cond_values_A, element_onset_ms(pos), time_win, bin_edges, n_time_bins, ...
            baseline_mu_fr, baseline_std_fr);
        trials_B(neuron_i).(position_classes{i}) = extract_bin_features( ...
            all_sdf, cond_value_all, cond_values_B, element_onset_ms(pos), time_win, bin_edges, n_time_bins, ...
            baseline_mu_fr, baseline_std_fr);
    end
end

%% ---- Step 2: cross-group generalization decode, per area ----

decoding_problems = {
    'Position_Auditory', auditory_neuron_idx;
    'Position_Frontal',  frontal_neuron_idx;
    };

crossgroup_results = struct();

for prob_i = 1:size(decoding_problems,1)

    label      = decoding_problems{prob_i,1};
    neuron_idx = decoding_problems{prob_i,2};

    fprintf('\n=== %s (cross-sequence-group generalization) ===\n', strrep(label,'_',' '));

    [boot_acc, boot_AtoB, boot_BtoA, perm_acc, n_used] = crossgroup_lda_decode( ...
        trials_A, trials_B, position_classes, neuron_idx, min_trials, n_pseudo_trial, n_time_bins, n_pcs, n_boot, n_perm);

    crossgroup_results.(label).boot_acc     = boot_acc;
    crossgroup_results.(label).boot_AtoB    = boot_AtoB;
    crossgroup_results.(label).boot_BtoA    = boot_BtoA;
    crossgroup_results.(label).perm_acc     = perm_acc;
    crossgroup_results.(label).neurons_used = n_used;

    obs_acc = median(boot_acc,'omitnan');
    p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));

    fprintf('Neurons usable (>= %d real trials/class, both group-sets): %d / %d\n', min_trials, n_used, length(neuron_idx));
    fprintf('Train A -> Test B: median = %.3f | Train B -> Test A: median = %.3f\n', ...
        median(boot_AtoB,'omitnan'), median(boot_BtoA,'omitnan'));
    fprintf('Combined (averaged) cross-group accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
    fprintf('Theoretical chance = 0.333 | Permutation null: median = %.3f, 95th pctile = %.3f\n', ...
        median(perm_acc,'omitnan'), prctile(perm_acc,95));
    fprintf('Permutation p-value (observed vs null): p = %.4f\n', p_perm);

end

%% ---- Step 3: plot ----

problem_names = fieldnames(crossgroup_results);
figure('Renderer','painters','Position',[100 100 500 450]); hold on

bar_x      = 1:length(problem_names);
bar_median = nan(size(bar_x));
err_lo     = nan(size(bar_x));
err_hi     = nan(size(bar_x));

for i = 1:length(problem_names)
    acc = crossgroup_results.(problem_names{i}).boot_acc;
    bar_median(i) = median(acc,'omitnan');
    err_lo(i)     = bar_median(i) - prctile(acc,2.5);
    err_hi(i)     = prctile(acc,97.5) - bar_median(i);
end

bar(bar_x, bar_median, 0.5, 'FaceColor',[0.4 0.6 0.8]);
errorbar(bar_x, bar_median, err_lo, err_hi, 'k', 'LineStyle','none','LineWidth',1.2)
plot(bar_x, 0.333*ones(size(bar_x)), 'rd','MarkerFaceColor','r')
set(gca,'XTick',bar_x,'XTickLabel',strrep(problem_names,'_',' '))
ylabel('Cross-sequence-group generalization accuracy')
ylim([0 1]); box off
title('Position decoding, trained on one pair of sequences / tested on the other')
legend({'Median accuracy','95% bootstrap CI','Theoretical chance'},'Location','best')

%% ---- Step 4: frontal vs auditory ----

fprintf('\nPosition decoding (cross-group generalization): frontal vs auditory\n');
bootstrap_compare(crossgroup_results.Position_Frontal.boot_acc, crossgroup_results.Position_Auditory.boot_acc);


%% ================= Local functions =================

function cond_values = lookup_cond_values_groupset(seq_letters, seq_cond_values, position, letter, allowed_groups)
% Same idea as lookup_cond_values in pca_lda_id_position_wholeseq_rigorous.m,
% but restricted to sequence groups in allowed_groups only -- this is what
% keeps Set A and Set B built from disjoint physical sequences.
cond_values = [];
for g = allowed_groups
    if strcmp(seq_letters{g,position}, letter)
        cond_values = [cond_values, seq_cond_values{g}]; %#ok<AGROW>
    end
end
end


function trial_bins = extract_bin_features(all_sdf, cond_value_all, cond_values, onset_ms, time_win, ...
    bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr)
% Identical to pca_lda_id_position_wholeseq_rigorous.m's version.

trial_idx = ismember(cond_value_all, cond_values);
if ~any(trial_idx)
    trial_bins = [];
    return
end

cols = 1000 + onset_ms + time_win;
trial_window = all_sdf(trial_idx, cols);

raw_bins = nan(size(trial_window,1), n_time_bins);
for bin_i = 1:n_time_bins
    raw_bins(:,bin_i) = nanmean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2);
end

valid_rows = all(~isnan(raw_bins),2);
raw_bins = raw_bins(valid_rows,:);

if isempty(raw_bins)
    trial_bins = [];
else
    trial_bins = (raw_bins - baseline_mu_fr) ./ baseline_std_fr;
end

end


function [boot_acc, boot_AtoB, boot_BtoA, perm_acc, n_used] = crossgroup_lda_decode( ...
    trials_A, trials_B, classes, neuron_idx, min_trials, n_pseudo_trial, n_time_bins, n_pcs, n_boot, n_perm)
% Trains on pseudo-populations built from Set A, tests on pseudo-
% populations built from Set B, and vice versa -- no within-set
% cross-validation is needed since the two sets are, by construction,
% disjoint physical trials from different sequences.

n_classes = length(classes);

usable       = true(length(neuron_idx),1);
trial_pool_A = cell(length(neuron_idx), n_classes);
trial_pool_B = cell(length(neuron_idx), n_classes);

for ni = 1:length(neuron_idx)
    for ci = 1:n_classes
        vA = trials_A(neuron_idx(ni)).(classes{ci});
        vB = trials_B(neuron_idx(ni)).(classes{ci});
        trial_pool_A{ni,ci} = vA;
        trial_pool_B{ni,ci} = vB;
        if isempty(vA) || size(vA,1) < min_trials || isempty(vB) || size(vB,1) < min_trials
            usable(ni) = false;
        end
    end
end

use_idx      = find(usable);
n_used       = length(use_idx);
trial_pool_A = trial_pool_A(use_idx,:);
trial_pool_B = trial_pool_B(use_idx,:);

boot_AtoB = nan(n_boot,1);
boot_BtoA = nan(n_boot,1);
for boot_i = 1:n_boot
    [X_A, y_A] = build_pseudo_population(trial_pool_A, n_pseudo_trial, n_time_bins);
    [X_B, y_B] = build_pseudo_population(trial_pool_B, n_pseudo_trial, n_time_bins);
    boot_AtoB(boot_i) = fit_test_accuracy(X_A, y_A, X_B, y_B, n_pcs);
    boot_BtoA(boot_i) = fit_test_accuracy(X_B, y_B, X_A, y_A, n_pcs);
end
boot_acc = mean([boot_AtoB, boot_BtoA], 2);

perm_AtoB = nan(n_perm,1);
perm_BtoA = nan(n_perm,1);
for perm_i = 1:n_perm
    [X_A, y_A] = build_pseudo_population(trial_pool_A, n_pseudo_trial, n_time_bins);
    [X_B, y_B] = build_pseudo_population(trial_pool_B, n_pseudo_trial, n_time_bins);
    y_A_shuf = y_A(randperm(length(y_A)));
    y_B_shuf = y_B(randperm(length(y_B)));
    perm_AtoB(perm_i) = fit_test_accuracy(X_A, y_A_shuf, X_B, y_B, n_pcs);
    perm_BtoA(perm_i) = fit_test_accuracy(X_B, y_B_shuf, X_A, y_A, n_pcs);
end
perm_acc = mean([perm_AtoB, perm_BtoA], 2);

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


function acc = fit_test_accuracy(X_train, y_train, X_test, y_test, n_pcs)

mu = mean(X_train,1);
[coeff, scores_train] = pca(X_train);
n_keep = min(n_pcs, size(coeff,2));
scores_train = scores_train(:,1:n_keep);
scores_test  = (X_test - mu) * coeff(:,1:n_keep);

lda_model = fitcdiscr(scores_train, y_train);
preds = predict(lda_model, scores_test);
acc = mean(preds == y_test);

end
