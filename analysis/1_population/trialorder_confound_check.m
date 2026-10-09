%% Trial-order confound check: does elapsed time in the session predict sequence identity?
%
% The manuscript states each testing trial's sequence was "randomly
% chosen" (Methods, Behavioural paradigm). If that randomization was
% genuinely trial-independent, knowing WHEN in the session a trial
% happened should carry no information about WHICH of the 4 grammatical
% sequences was played. If it does, that's a plausible mundane explanation
% for the confound found in sequence_confound_check.m (frontal 89% /
% auditory 67% decoding of sequence identity from a physically identical
% "A" stimulus) -- e.g. a pseudorandomization scheme that avoids immediate
% repeats, or blocking, could let simple session-time drift in neural
% state masquerade as anticipatory/predictive coding.
%
% This check uses NO neural data at all -- only each session's event
% table (trial_n, cond_value, cond_label) -- so it is fast and, if trial
% order turns out NOT to predict sequence identity, it directly rules out
% the mundane explanation rather than just failing to find it.
%
% Assumes: dirs, spike_log already in the workspace (spike_log.session is
% used only to get the list of unique recording sessions).

clear pooled_trial_n pooled_group

rng(20,'twister')

n_boot     = 300;
n_perm     = 300;
k_folds    = 5;
cv_repeats = 4;

seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};   % same 4 sequence groups as sequence_confound_check.m

unique_sessions = unique(spike_log.session);
n_sessions = length(unique_sessions);

pooled_trial_n = [];
pooled_group   = [];
pooled_session = {};

for session_i = 1:n_sessions

    try
        event_table_in = load(fullfile(dirs.mat_data, [unique_sessions{session_i} '.mat']), 'event_table');
    catch
        continue
    end

    event_table = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);

    trial_n_session    = event_table.trial_n(nonviol_mask);
    cond_value_session = event_table.cond_value(nonviol_mask);

    group_session = nan(size(cond_value_session));
    for g = 1:4
        group_session(ismember(cond_value_session, seq_cond_values{g})) = g;
    end

    valid = ~isnan(group_session) & ~isnan(trial_n_session);
    trial_n_session = trial_n_session(valid);
    group_session   = group_session(valid);

    if length(trial_n_session) < 8   % need at least a couple of trials per group to be meaningful
        continue
    end

    % z-score trial number WITHIN session so sessions of different length
    % and different absolute trial-count ranges don't dominate the pooled fit
    trial_n_z = (trial_n_session - mean(trial_n_session)) ./ std(trial_n_session);

    pooled_trial_n = [pooled_trial_n; trial_n_z]; %#ok<AGROW>
    pooled_group   = [pooled_group; group_session]; %#ok<AGROW>
    pooled_session = [pooled_session; repmat(unique_sessions(session_i), length(trial_n_z), 1)]; %#ok<AGROW>

end

fprintf('Pooled nonviolation trials with a valid sequence group, across %d sessions: %d\n', ...
    length(unique(pooled_session)), length(pooled_group));

%% Decode sequence group (4-way, chance = 0.25) from within-session trial order alone

boot_acc = nan(n_boot,1);
for boot_i = 1:n_boot
    idx = randi(length(pooled_group), length(pooled_group), 1);   % resample trials with replacement
    boot_acc(boot_i) = nested_cv_accuracy_1d(pooled_trial_n(idx), pooled_group(idx), k_folds, cv_repeats);
end

perm_acc = nan(n_perm,1);
for perm_i = 1:n_perm
    idx = randi(length(pooled_group), length(pooled_group), 1);
    y_perm = pooled_group(idx);
    y_perm = y_perm(randperm(length(y_perm)));   % break the trial-order/group relationship
    perm_acc(perm_i) = nested_cv_accuracy_1d(pooled_trial_n(idx), y_perm, k_folds, cv_repeats);
end

obs_acc = median(boot_acc,'omitnan');
p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));

fprintf('\n=== Trial order -> sequence identity (no neural data, behavioural log only) ===\n');
fprintf('Nested-CV accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
    obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
fprintf('Theoretical chance = 0.250 | Permutation null: median = %.3f, 95th pctile = %.3f\n', ...
    median(perm_acc,'omitnan'), prctile(perm_acc,95));
fprintf('Permutation p-value: p = %.4f\n', p_perm);

if obs_acc <= prctile(perm_acc,95)
    fprintf(['--> Trial order does NOT predict sequence identity above the permutation null.\n' ...
             '    Randomization looks intact: the frontal identical-stimulus confound is NOT explained\n' ...
             '    by simple within-session trial-order drift, and is a more genuinely interesting result.\n']);
else
    fprintf(['--> Trial order DOES predict sequence identity above the permutation null.\n' ...
             '    Sequence assignment was not fully independent of trial order -- this is a plausible\n' ...
             '    (at least partial) mundane explanation for the neural confound found earlier.\n']);
end

%% ================= Local function =================

function acc = nested_cv_accuracy_1d(x, y, k_folds, cv_repeats)
% Same repeated stratified K-fold CV logic as the other rigorous scripts,
% but with a single scalar predictor (trial order) instead of a PCA-
% reduced neural population vector -- no PCA step needed for one feature.

fold_acc = nan(k_folds*cv_repeats,1);
f = 1;
for rep_i = 1:cv_repeats
    cvp = cvpartition(y, 'KFold', k_folds);
    for fold_i = 1:k_folds
        train_idx = training(cvp, fold_i);
        test_idx  = test(cvp, fold_i);

        lda_model = fitcdiscr(x(train_idx), y(train_idx));
        preds = predict(lda_model, x(test_idx));
        fold_acc(f) = mean(preds == y(test_idx));
        f = f + 1;
    end
end
acc = mean(fold_acc,'omitnan');

end
