%% Single-neuron and PCA-free population decoding of ordinal position, with and without the sequence-group confound
%
% Purpose: a fair, pre-specified test of whether position information in
% FIRING RATES (no PCA anywhere) is greater in frontal than auditory cortex,
% at the level of (1) every individual neuron and (2) the pseudo-population.
%
% Why two decoders per neuron: the design has only 4 grammatical sequences,
% so "which position is C at" is entangled with "which sequence is this" (and
% with trial-N-1 carryover, which is decodable at 0.985 in frontal). A decoder
% evaluated by ordinary cross-validation within the pooled trials can use that
% entanglement. So every neuron is decoded two ways, with the identical
% classifier and features:
%
%   WITHIN-GROUP CV  -- 5-fold CV over trials pooled across all four
%                       sequences. Position labels are aliased with sequence
%                       group here, so this is the "what the raw data shows"
%                       number and includes any sequence-context signal.
%   CROSS-GROUP      -- train on trials from sequences {1,2}, test on trials
%                       from sequences {3,4} (and vice versa). Train and test
%                       come from physically different sequences, so only a
%                       position code that GENERALISES ACROSS SEQUENCES can
%                       score above chance. This is exactly the property the
%                       Supplementary Discussion claims ("ordinal decoding
%                       generalized across multiple well-learned sequences
%                       whose elements differed at matched ordinal
%                       positions"), so it is a test of the manuscript's own
%                       stated claim, not a stricter standard imposed on it.
%
% The gap between the two is itself informative: if WITHIN is high in frontal
% while CROSS is at chance, the apparent effect in the raw data is
% sequence-context information, not ordinal position.
%
% PRE-SPECIFIED READING (written before any result exists):
%   The position claim is supported only if, on the CROSS-GROUP measure,
%     (i)  the fraction of significant neurons is higher in frontal than
%          auditory (Fisher exact, p < 0.05), AND
%     (ii) the frontal fraction exceeds the 5% false-positive base rate
%          (binomial, p < 0.05), AND
%     (iii) the PCA-free population cross-group accuracy in frontal exceeds
%          its label-permutation null (p < 0.05) and exceeds auditory
%          (bootstrap p < 0.05).
%   WITHIN-GROUP results alone do not count as support, for the reasons above.
%
% Features: single-trial, sequence-onset-anchored (the same data source behind
% the Fig 2A/B trajectory result, so any slow ramp that carries elapsed-time
% information is available to the decoder), 4 time bins over 0-413 ms after
% each element onset, z-scored to a global pre-sequence baseline. Positions
% 2, 3, 5 of identity C (position 4 is produced by only one sequence and
% cannot be separated from sequence identity at all).
%
% Classifier: diagonal-covariance nearest-centroid (equivalent to diagonal
% LDA with equal priors), applied directly to the raw z-scored features, no
% dimensionality reduction. Scored with balanced accuracy (chance = 1/3).
%
% Assumes in the workspace (as set up by aglt_analysis_main.m): spike_log,
% dirs, auditory_neuron_idx, frontal_neuron_idx. Also uses bootstrap_compare.m.
%
% Runtime: loading all neurons' sequence-onset SDFs takes a few minutes; the
% single-neuron permutation tests take a few more; the population bootstrap
% and permutation loops are the longest part. Reduce n_boot/n_perm for a
% quick first pass.

clear trials_A trials_B

rng(22,'twister')

%% ---- User-set parameters ----
element_onset_ms = [0 563 1126 1688 2252];
time_win         = 0:413;
n_time_bins      = 4;

min_trials_sn    = 5;      % per class, per set, per neuron (single-neuron decoder)
min_trials_pop   = 3;      % per class, per set, per neuron (pseudo-population)
n_perm_sn        = 200;    % label permutations per neuron
k_folds_sn       = 5;      % within-group CV folds
cv_repeats_sn    = 2;      % repeats of the within-group fold assignment
alpha_sn         = 0.05;   % per-neuron significance threshold

n_pseudo_trial   = 15;     % pseudo-trials per class per pseudo-population draw
n_boot           = 300;    % population bootstrap iterations
n_perm           = 300;    % population label-permutation iterations
n_neuron_boot    = 1000;   % neuron-resampling bootstrap for single-neuron medians

bin_edges = round(linspace(0, length(time_win), n_time_bins+1));

seq_letters = { ...
    'A','C','G','F','C'; ...   % cond_value 1 & 5   (Set A)
    'A','D','C','G','F'; ...   % cond_value 2 & 6   (Set A)
    'A','C','F','C','G'; ...   % cond_value 3 & 7   (Set B)
    'A','D','C','F','C'};      % cond_value 4 & 8   (Set B)
seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};

groupset_A = [1 2];
groupset_B = [3 4];

position_list    = [2 3 5];
position_classes = {'position_2','position_3','position_5'};
K = numel(position_classes);

fprintf('Cross-group assignment (position -> sequence used in Set A / Set B):\n');
for i = 1:length(position_list)
    pos  = position_list(i);
    cand = find(strcmp(seq_letters(:,pos), 'C'))';
    gA = intersect(cand, groupset_A);
    gB = intersect(cand, groupset_B);
    fprintf('  position %d : Set A = seq%d, Set B = seq%d\n', pos, gA(1), gB(1));
end

%% ---- Step 1: single-trial, sequence-onset-anchored, baseline z-scored bin features ----

n_neurons = size(spike_log,1);
trials_A(n_neurons,1) = struct();
trials_B(n_neurons,1) = struct();

for neuron_i = 1:n_neurons

    if mod(neuron_i,100) == 0
        fprintf('Loading neuron %i of %i\n', neuron_i, n_neurons);
    end

    try
        sdf_in = load(fullfile(dirs.root, 'data', 'spike', ...
            [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
        event_table_in = load(fullfile(dirs.mat_data, ...
            [spike_log.session{neuron_i} '.mat']), 'event_table');
    catch
        continue
    end

    event_table  = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);

    baseline_window = 800:1000;
    all_sdf         = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    cond_value_all  = event_table.cond_value(nonviol_mask);

    baseline_mu_fr  = mean(mean(all_sdf(:,baseline_window),2,'omitnan'),'omitnan');
    baseline_std_fr = std(mean(all_sdf(:,baseline_window),2,'omitnan'),'omitnan');

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

%% ---- Step 2: single-neuron decoding (cross-group and within-group), every neuron ----

sn_cross_acc  = nan(n_neurons,1);
sn_cross_p    = nan(n_neurons,1);
sn_within_acc = nan(n_neurons,1);
sn_within_p   = nan(n_neurons,1);

for neuron_i = 1:n_neurons

    if mod(neuron_i,200) == 0
        fprintf('Single-neuron decoding: neuron %i of %i\n', neuron_i, n_neurons);
    end

    [X_A, y_A, okA] = stack_classes(trials_A(neuron_i), position_classes, min_trials_sn);
    [X_B, y_B, okB] = stack_classes(trials_B(neuron_i), position_classes, min_trials_sn);
    if ~(okA && okB)
        continue
    end

    % Cross-group (A->B and B->A)
    obs_c = crossgroup_acc(X_A, y_A, X_B, y_B, K, y_A, y_B);
    null_c = nan(n_perm_sn,1);
    for perm_i = 1:n_perm_sn
        null_c(perm_i) = crossgroup_acc(X_A, y_A, X_B, y_B, K, ...
            y_A(randperm(numel(y_A))), y_B(randperm(numel(y_B))));
    end
    sn_cross_acc(neuron_i) = obs_c;
    sn_cross_p(neuron_i)   = (1 + sum(null_c >= obs_c)) / (1 + n_perm_sn);

    % Within-group CV over the pooled trials (sequence group NOT controlled)
    X_all = [X_A; X_B];
    y_all = [y_A; y_B];
    obs_w = within_cv_acc(X_all, y_all, K, k_folds_sn, cv_repeats_sn);
    null_w = nan(n_perm_sn,1);
    for perm_i = 1:n_perm_sn
        null_w(perm_i) = within_cv_acc(X_all, y_all(randperm(numel(y_all))), K, k_folds_sn, cv_repeats_sn);
    end
    sn_within_acc(neuron_i) = obs_w;
    sn_within_p(neuron_i)   = (1 + sum(null_w >= obs_w)) / (1 + n_perm_sn);
end

%% ---- Step 3: summarise the single-neuron results by area ----

area_names = {'Auditory','Frontal'};
area_idx   = {auditory_neuron_idx, frontal_neuron_idx};

n_used   = nan(1,2);
k_cross  = nan(1,2);
k_within = nan(1,2);
acc_cross_by_area  = cell(1,2);
acc_within_by_area = cell(1,2);

fprintf('\n=== Single-neuron position decoding (positions 2/3/5 of C; chance = 0.333) ===\n');
for a = 1:2
    idx = area_idx{a};
    ok  = idx(~isnan(sn_cross_acc(idx)));
    n_used(a)   = numel(ok);
    k_cross(a)  = sum(sn_cross_p(ok)  < alpha_sn);
    k_within(a) = sum(sn_within_p(ok) < alpha_sn);
    acc_cross_by_area{a}  = sn_cross_acc(ok);
    acc_within_by_area{a} = sn_within_acc(ok);

    p_bin_cross  = 1 - binocdf(k_cross(a)  - 1, n_used(a), alpha_sn);
    p_bin_within = 1 - binocdf(k_within(a) - 1, n_used(a), alpha_sn);

    fprintf('\n%s: %d / %d neurons usable (>= %d trials per class in both sequence sets)\n', ...
        area_names{a}, n_used(a), numel(idx), min_trials_sn);
    fprintf('  CROSS-group : median acc = %.3f | %d neurons p<%.2f (%.1f%%, base rate %.0f%%), binomial p = %.4f\n', ...
        median(acc_cross_by_area{a}), k_cross(a), alpha_sn, 100*k_cross(a)/n_used(a), 100*alpha_sn, p_bin_cross);
    fprintf('  WITHIN-group: median acc = %.3f | %d neurons p<%.2f (%.1f%%, base rate %.0f%%), binomial p = %.4f\n', ...
        median(acc_within_by_area{a}), k_within(a), alpha_sn, 100*k_within(a)/n_used(a), 100*alpha_sn, p_bin_within);
end

fprintf('\n--- Frontal vs auditory, single neurons ---\n');

[~, p_fisher_cross]  = fishertest([k_cross(2)  n_used(2)-k_cross(2);  k_cross(1)  n_used(1)-k_cross(1)]);
[~, p_fisher_within] = fishertest([k_within(2) n_used(2)-k_within(2); k_within(1) n_used(1)-k_within(1)]);
fprintf('PRIMARY  CROSS-group  fraction significant: frontal %.1f%% vs auditory %.1f%%, Fisher exact p = %.4f\n', ...
    100*k_cross(2)/n_used(2), 100*k_cross(1)/n_used(1), p_fisher_cross);
fprintf('context  WITHIN-group fraction significant: frontal %.1f%% vs auditory %.1f%%, Fisher exact p = %.4f\n', ...
    100*k_within(2)/n_used(2), 100*k_within(1)/n_used(1), p_fisher_within);

rng(23,'twister')
for measure_i = 1:2
    if measure_i == 1
        accF = acc_cross_by_area{2};  accA = acc_cross_by_area{1};  mname = 'CROSS-group';
    else
        accF = acc_within_by_area{2}; accA = acc_within_by_area{1}; mname = 'WITHIN-group';
    end
    nF = numel(accF); nA = numel(accA);
    mF = nan(n_neuron_boot,1); mA = nan(n_neuron_boot,1);
    for b = 1:n_neuron_boot
        mF(b) = median(accF(randi(nF,nF,1)));
        mA(b) = median(accA(randi(nA,nA,1)));
    end
    fprintf('\n%s per-neuron accuracy, frontal (median %.3f) vs auditory (median %.3f): rank-sum p = %.4f\n', ...
        mname, median(accF), median(accA), ranksum(accF, accA));
    fprintf('Neuron-resampling bootstrap median difference (frontal - auditory):\n');
    bootstrap_compare(mF, mA);
end

fprintf('\nWITHIN minus CROSS median accuracy (inflation attributable to sequence context): auditory %.3f, frontal %.3f\n', ...
    median(acc_within_by_area{1}) - median(acc_cross_by_area{1}), ...
    median(acc_within_by_area{2}) - median(acc_cross_by_area{2}));

%% ---- Step 4: PCA-free population decoding, cross-group (pseudo-population) ----

pop_results = struct();
for a = 1:2
    fprintf('\n=== Position, %s, PCA-free population, cross-sequence-group ===\n', area_names{a});
    [boot_acc, perm_acc, n_pop] = crossgroup_pop_decode(trials_A, trials_B, position_classes, area_idx{a}, ...
        min_trials_pop, n_pseudo_trial, n_time_bins, n_boot, n_perm);
    pop_results.(area_names{a}).boot_acc = boot_acc;
    pop_results.(area_names{a}).perm_acc = perm_acc;
    obs_acc = median(boot_acc,'omitnan');
    p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));
    fprintf('Neurons usable (>= %d trials/class in both sets): %d / %d\n', min_trials_pop, n_pop, numel(area_idx{a}));
    fprintf('Cross-group accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
    fprintf('Chance = 0.333 | Permutation null: median = %.3f, 95th pctile = %.3f | permutation p = %.4f\n', ...
        median(perm_acc,'omitnan'), prctile(perm_acc,95), p_perm);
end

fprintf('\nPRIMARY  Population cross-group position decoding: frontal vs auditory\n');
bootstrap_compare(pop_results.Frontal.boot_acc, pop_results.Auditory.boot_acc);

%% ---- Step 5: plots (standard MATLAB plotting) ----

figure('Renderer','painters','Position',[100 100 650 400]); hold on
group_vals   = {acc_within_by_area{1}, acc_within_by_area{2}, acc_cross_by_area{1}, acc_cross_by_area{2}};
group_names  = {'Aud within','Fro within','Aud cross','Fro cross'};
group_colors = [0.5 0.6 0.8; 0.9 0.6 0.5; 0.2 0.4 0.7; 0.8 0.3 0.2];
for g = 1:4
    v = group_vals{g};
    jitter_x = g + 0.15*(rand(size(v))-0.5);
    scatter(jitter_x, v, 8, group_colors(g,:), 'filled', 'MarkerFaceAlpha', 0.3);
    plot(g + [-0.25 0.25], median(v)*[1 1], 'k-', 'LineWidth', 2);
end
plot([0.5 4.5], [1/3 1/3], 'r--');
set(gca,'XTick',1:4,'XTickLabel',group_names); ylabel('Single-neuron balanced accuracy'); box off
title('Single-neuron position decoding: within-group CV vs cross-sequence-group (black = median)')

figure('Renderer','painters','Position',[100 100 450 400]); hold on
pop_names = {'Auditory','Frontal'};
for a = 1:2
    acc = pop_results.(pop_names{a}).boot_acc;
    m = median(acc,'omitnan');
    bar(a, m, 0.5, 'FaceColor', group_colors(a+2,:));
    errorbar(a, m, m - prctile(acc,2.5), prctile(acc,97.5) - m, 'k', 'LineWidth', 1.2);
    plot(a + [-0.3 0.3], prctile(pop_results.(pop_names{a}).perm_acc,95)*[1 1], 'k--');
end
plot([0.5 2.5], [1/3 1/3], 'r--');
set(gca,'XTick',1:2,'XTickLabel',pop_names); ylabel('Cross-group accuracy'); ylim([0 1]); box off
title('PCA-free population cross-group decoding (dashed black = permutation 95th pct)')


%% ================= Local functions =================

function cond_values = lookup_cond_values_groupset(seq_letters, seq_cond_values, position, letter, allowed_groups)
cond_values = [];
for g = allowed_groups
    if strcmp(seq_letters{g,position}, letter)
        cond_values = [cond_values, seq_cond_values{g}]; %#ok<AGROW>
    end
end
end


function trial_bins = extract_bin_features(all_sdf, cond_value_all, cond_values, onset_ms, time_win, ...
    bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr)

trial_idx = ismember(cond_value_all, cond_values);
if ~any(trial_idx)
    trial_bins = [];
    return
end

cols = 1000 + onset_ms + time_win;
trial_window = all_sdf(trial_idx, cols);

raw_bins = nan(size(trial_window,1), n_time_bins);
for bin_i = 1:n_time_bins
    raw_bins(:,bin_i) = mean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2, 'omitnan');
end

valid_rows = all(~isnan(raw_bins),2);
raw_bins = raw_bins(valid_rows,:);

if isempty(raw_bins)
    trial_bins = [];
else
    trial_bins = (raw_bins - baseline_mu_fr) ./ baseline_std_fr;
end

end


function [X, y, ok] = stack_classes(s, classes, min_trials)
% Stack one neuron's per-class trial matrices into X (trials x bins) and a
% class-index vector y. ok is false if any class is missing or too sparse.
X = []; y = []; ok = true;
for ci = 1:numel(classes)
    if ~isfield(s, classes{ci}) || isempty(s.(classes{ci})) || size(s.(classes{ci}),1) < min_trials
        ok = false;
        return
    end
    v = s.(classes{ci});
    X = [X; v];                              %#ok<AGROW>
    y = [y; ci*ones(size(v,1),1)];           %#ok<AGROW>
end
end


function pred = diag_nc_predict(Xtr, ytr, Xte, K)
% Diagonal-covariance nearest-centroid classifier (diagonal LDA, equal
% priors) on raw features. No dimensionality reduction.
mu = nan(K, size(Xtr,2));
for k = 1:K
    mu(k,:) = mean(Xtr(ytr==k,:), 1);
end
resid = Xtr - mu(ytr,:);
v = sum(resid.^2, 1) ./ max(size(Xtr,1) - K, 1);
v = max(v, 1e-3);
d = nan(size(Xte,1), K);
for k = 1:K
    d(:,k) = sum(((Xte - mu(k,:)).^2) ./ v, 2);
end
[~, pred] = min(d, [], 2);
end


function ba = balanced_acc(pred, y, K)
r = nan(1,K);
for k = 1:K
    m = (y == k);
    r(k) = mean(pred(m) == k);
end
ba = mean(r, 'omitnan');
end


function acc = crossgroup_acc(X_A, y_A, X_B, y_B, K, yA_train, yB_train)
% Train on A (labels yA_train) -> test on B (true labels y_B), and vice
% versa; return the mean of the two balanced accuracies. For the permutation
% null, pass shuffled yA_train / yB_train and the TRUE y_A / y_B.
pAB = diag_nc_predict(X_A, yA_train, X_B, K);
pBA = diag_nc_predict(X_B, yB_train, X_A, K);
acc = mean([balanced_acc(pAB, y_B, K), balanced_acc(pBA, y_A, K)]);
end


function acc = within_cv_acc(X, y, K, k_folds, repeats)
% Repeated stratified k-fold CV over the pooled trials (sequence group is
% NOT controlled here -- this is the confound-inclusive comparison).
accs = nan(repeats,1);
for r = 1:repeats
    fold = stratified_folds(y, k_folds);
    pred = nan(size(y));
    for f = 1:k_folds
        te = (fold == f);
        pred(te) = diag_nc_predict(X(~te,:), y(~te), X(te,:), K);
    end
    accs(r) = balanced_acc(pred, y, K);
end
acc = mean(accs);
end


function fold = stratified_folds(y, k)
fold = zeros(size(y));
for c = unique(y)'
    idx = find(y == c);
    idx = idx(randperm(numel(idx)));
    fold(idx) = mod(0:numel(idx)-1, k)' + 1;
end
end


function [boot_acc, perm_acc, n_used] = crossgroup_pop_decode(trials_A, trials_B, classes, neuron_idx, ...
    min_trials, n_pseudo_trial, n_time_bins, n_boot, n_perm)
% PCA-free pseudo-population cross-group decode: train on pseudo-populations
% built from Set A, test on Set B, and vice versa.

n_classes = length(classes);
K = n_classes;

usable       = true(length(neuron_idx),1);
trial_pool_A = cell(length(neuron_idx), n_classes);
trial_pool_B = cell(length(neuron_idx), n_classes);

for ni = 1:length(neuron_idx)
    for ci = 1:n_classes
        if ~isfield(trials_A(neuron_idx(ni)), classes{ci}) || ~isfield(trials_B(neuron_idx(ni)), classes{ci})
            usable(ni) = false;
            continue
        end
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
    boot_AtoB(boot_i) = balanced_acc(diag_nc_predict(X_A, y_A, X_B, K), y_B, K);
    boot_BtoA(boot_i) = balanced_acc(diag_nc_predict(X_B, y_B, X_A, K), y_A, K);
end
boot_acc = mean([boot_AtoB, boot_BtoA], 2);

perm_AtoB = nan(n_perm,1);
perm_BtoA = nan(n_perm,1);
for perm_i = 1:n_perm
    [X_A, y_A] = build_pseudo_population(trial_pool_A, n_pseudo_trial, n_time_bins);
    [X_B, y_B] = build_pseudo_population(trial_pool_B, n_pseudo_trial, n_time_bins);
    perm_AtoB(perm_i) = balanced_acc(diag_nc_predict(X_A, y_A(randperm(length(y_A))), X_B, K), y_B, K);
    perm_BtoA(perm_i) = balanced_acc(diag_nc_predict(X_B, y_B(randperm(length(y_B))), X_A, K), y_A, K);
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
