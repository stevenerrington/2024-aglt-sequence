%% Position decoding with IDENTITY and PREDECESSOR held constant, tested across
%% sequences that never overlap (grammatical vs. violation sequences)
%
% WHY THIS IS A DIFFERENT KIND OF TEST
% Every earlier analysis used only the 4 grammatical sequences, where "which
% position is this element at" is aliased with "which sequence is this" (and
% with trial N-1 carryover). The 8 violation sequences (cond_value 9-16) contain
% the SAME elements in new orders. That gives minimal pairs: the same letter,
% preceded by the same letter, at two different ordinal positions, in several
% different sequences, so position can be separated from sequence identity by
% design instead of by a statistical correction.
%
%   Primary context  : C preceded by F, position 4 vs position 5
%       position 4, grammatical-equivalent : cond 3,7 (A C F C G) and cond 15
%                    (A C F C F, identical up to the violation at position 5)
%       position 5, grammatical-equivalent : cond 1,5 (A C G F C), 4,8 (A D C F C)
%       position 4, post-violation          : cond 9 (A G F C D), 13 (A D F C G)
%       position 5, post-violation          : cond 10 (A C D F C), 11 (A G C F C),
%                                             14 (A D G F C)
%   Secondary context: G preceded by C, position 4 vs position 5 (thin: one
%       post-violation sequence per position; reported as replication only).
%
% THE TEST. Train a position classifier (pos_lo vs pos_hi) on one set
% (grammatical-equivalent) and test on the other (post-violation), and vice
% versa. The two sets share no sequence, so sequence-identity and trial-N-1
% carryover cannot transfer. Identity and immediate predecessor are identical
% in both classes by construction. Chance = 0.5.
%
% PRE-SPECIFIED READING (written before any result exists). The ordinal-position
% claim is supported only if, for the PRIMARY context on the CROSS-SET measure,
%   (i)   frontal PCA-free population accuracy exceeds its permutation 95th
%         percentile, AND
%   (ii)  frontal population accuracy exceeds auditory (bootstrap two-sided
%         p < 0.05, frontal higher), AND
%   (iii) the fraction of frontal neurons significant on the cross-set measure
%         exceeds the 5% base rate (binomial p < 0.05) AND exceeds auditory
%         (Fisher exact p < 0.05, frontal higher).
% The secondary context cannot rescue a failure of the primary. Within-set
% results (which keep position aliased with sequence) do not count as support.
%
% LIMITS THAT NO ANALYSIS OF THIS DATASET CAN REMOVE
%   - Ordinal position is perfectly confounded with elapsed time since sequence
%     onset (and time to reward). A positive result means "a position/elapsed-
%     time signal that generalises across sequences", not position per se.
%   - In the post-violation set, time since the violation covaries with position
%     (pos4 cells are 1-2 elements after the violation, pos5 cells 2-3). Because
%     transfer requires the same neural axis to separate classes in BOTH sets,
%     and grammatical-equivalent cells contain no violation, that confound can
%     only help if the surprise-recency axis coincides with the position axis.
%     The G-vs-V context-shift control below quantifies how large the
%     violation-related shift is.
%
% CONTROLS REPORTED ALONGSIDE
%   - Within-set split-half decoding (position still aliased with sequence):
%     shows whether the pipeline can detect information at all.
%   - Context shift: decode grammatical-equivalent vs post-violation at a FIXED
%     position, identity and predecessor (how large is the surprise shift).
%
% Classifier / features are identical to singleunit_position_crossgroup.m:
% baseline z-scored single-trial firing-rate bins (4 bins, 0-413 ms after
% element onset), diagonal nearest-centroid, balanced accuracy, no PCA.
%
% Assumes in the workspace (as set up by aglt_analysis_main.m): spike_log, dirs,
% auditory_neuron_idx, frontal_neuron_idx. Uses bootstrap_compare.m.
%
% Assumption to verify: cond_value 9-12 = Multiviolation_1-4 and 13-16 =
% Single_Rule_break_1-4, in the order of 'sequence description.xlsx'
% (matches aglt_save_spikes.m: 13,14 -> violation at 1127 ms, 15,16 -> 2253 ms).

clear cdat trG trV

rng(31,'twister')

%% ---- User-set parameters ----
element_onset_ms = [0 563 1126 1688 2252];
time_win         = 0:413;
n_time_bins      = 4;
bin_edges        = round(linspace(0, length(time_win), n_time_bins+1));

min_trials_sn    = 4;      % per cell (class x set), per neuron: single-neuron decoder
min_trials_pop   = 3;      % per cell, per neuron: pseudo-population
n_perm_sn        = 200;
k_folds_sn       = 4;
cv_repeats_sn    = 2;
alpha_sn         = 0.05;

n_pseudo_trial   = 15;
n_boot           = 300;
n_perm           = 300;
n_neuron_boot    = 1000;

% Letters of all 16 conditions (rows = cond_value 1..16)
cond_letters = { ...
    'A','C','G','F','C'; ...   %  1
    'A','D','C','G','F'; ...   %  2
    'A','C','F','C','G'; ...   %  3
    'A','D','C','F','C'; ...   %  4
    'A','C','G','F','C'; ...   %  5  (same sequence as 1)
    'A','D','C','G','F'; ...   %  6  (same as 2)
    'A','C','F','C','G'; ...   %  7  (same as 3)
    'A','D','C','F','C'; ...   %  8  (same as 4)
    'A','G','F','C','D'; ...   %  9  Multiviolation_1
    'A','C','D','F','C'; ...   % 10  Multiviolation_2
    'A','G','C','F','C'; ...   % 11  Multiviolation_3
    'A','F','C','G','C'; ...   % 12  Multiviolation_4
    'A','D','F','C','G'; ...   % 13  Single_Rule_break_1
    'A','D','G','F','C'; ...   % 14  Single_Rule_break_2
    'A','C','F','C','F'; ...   % 15  Single_Rule_break_3
    'A','D','C','G','C'};      % 16  Single_Rule_break_4

% First position at which each condition departs from the grammar (NaN = none).
% Elements at positions BEFORE this are identical in history to a grammatical sequence.
first_viol_pos = [nan(1,8), 2 3 2 2 3 3 5 5];

% Contexts: the pair is (letter | predecessor) at pos_lo vs pos_hi. First row = PRIMARY.
contexts = struct( ...
    'letter', {'C','G'}, ...
    'prev',   {'F','C'}, ...
    'pos_lo', {4, 4}, ...
    'pos_hi', {5, 5});
n_ctx = numel(contexts);

ctx_cond = struct('gLo',cell(1,n_ctx),'gHi',[],'vLo',[],'vHi',[]);
fprintf('Minimal-pair design (identity | predecessor held constant):\n');
for k = 1:n_ctx
    [ctx_cond(k).gLo, ctx_cond(k).vLo] = context_cells(cond_letters, first_viol_pos, ...
        contexts(k).letter, contexts(k).prev, contexts(k).pos_lo);
    [ctx_cond(k).gHi, ctx_cond(k).vHi] = context_cells(cond_letters, first_viol_pos, ...
        contexts(k).letter, contexts(k).prev, contexts(k).pos_hi);
    fprintf('  context %d: %s after %s, position %d vs %d%s\n', k, contexts(k).letter, contexts(k).prev, ...
        contexts(k).pos_lo, contexts(k).pos_hi, ternary_str(k==1,'  (PRIMARY)','  (secondary)'));
    fprintf('     pos %d: grammatical-equivalent conds %s | post-violation conds %s\n', contexts(k).pos_lo, ...
        mat2str(ctx_cond(k).gLo), mat2str(ctx_cond(k).vLo));
    fprintf('     pos %d: grammatical-equivalent conds %s | post-violation conds %s\n', contexts(k).pos_hi, ...
        mat2str(ctx_cond(k).gHi), mat2str(ctx_cond(k).vHi));
end

%% ---- Step 1: single-trial, sequence-onset-anchored, baseline z-scored bin features ----

n_neurons = size(spike_log,1);
cdat = cell(n_neurons, n_ctx);
viol_filter_mode = nan(n_neurons,1);   % 1 = rewarded violation trials only, 0 = all violation trials

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

    event_table = event_table_in.event_table;
    lab         = event_table.cond_label;
    rewarded    = ~isnan(event_table.rewardOnset_ms);

    is_nonviol  = strcmp(lab,'nonviol') & rewarded;
    is_viol_all = strcmp(lab,'viol');
    is_viol_rew = is_viol_all & rewarded;
    if any(is_viol_rew)
        is_viol = is_viol_rew;
        viol_filter_mode(neuron_i) = 1;
    elseif any(is_viol_all)
        is_viol = is_viol_all;
        viol_filter_mode(neuron_i) = 0;
    else
        is_viol = false(size(is_viol_all));
    end
    use_mask = is_nonviol | is_viol;

    baseline_window = 800:1000;
    all_sdf         = sdf_in.sdf.sequenceOnset(use_mask,:);
    cond_value_all  = event_table.cond_value(use_mask);

    baseline_mu_fr  = mean(mean(all_sdf(:,baseline_window),2,'omitnan'),'omitnan');
    baseline_std_fr = std(mean(all_sdf(:,baseline_window),2,'omitnan'),'omitnan');

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue
    end

    for k = 1:n_ctx
        d = struct();
        d.gLo = extract_bin_features(all_sdf, cond_value_all, ctx_cond(k).gLo, element_onset_ms(contexts(k).pos_lo), ...
            time_win, bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr);
        d.gHi = extract_bin_features(all_sdf, cond_value_all, ctx_cond(k).gHi, element_onset_ms(contexts(k).pos_hi), ...
            time_win, bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr);
        d.vLo = extract_bin_features(all_sdf, cond_value_all, ctx_cond(k).vLo, element_onset_ms(contexts(k).pos_lo), ...
            time_win, bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr);
        d.vHi = extract_bin_features(all_sdf, cond_value_all, ctx_cond(k).vHi, element_onset_ms(contexts(k).pos_hi), ...
            time_win, bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr);
        cdat{neuron_i,k} = d;
    end
end

fprintf('\nViolation trial filter: %d neurons used rewarded violation trials only, %d fell back to all violation trials, %d had none.\n', ...
    sum(viol_filter_mode==1), sum(viol_filter_mode==0), sum(isnan(viol_filter_mode)));

area_names = {'Auditory','Frontal'};
area_idx   = {auditory_neuron_idx, frontal_neuron_idx};

%% ---- Step 2: single-neuron decoding of position, cross-set and within-set ----

sn_cross_acc  = nan(n_neurons, n_ctx);
sn_cross_p    = nan(n_neurons, n_ctx);
sn_within_acc = nan(n_neurons, n_ctx);
sn_within_p   = nan(n_neurons, n_ctx);
sn_cell_n     = nan(n_neurons, 4, n_ctx);

for k = 1:n_ctx
    for neuron_i = 1:n_neurons

        if mod(neuron_i,200) == 0
            fprintf('Single-neuron decoding, context %d: neuron %i of %i\n', k, neuron_i, n_neurons);
        end

        d = cdat{neuron_i,k};
        if isempty(d); continue; end
        nn = [size(d.gLo,1) size(d.gHi,1) size(d.vLo,1) size(d.vHi,1)];
        sn_cell_n(neuron_i,:,k) = nn;
        if any(nn < min_trials_sn); continue; end

        [X_G, y_G] = stack2(d.gLo, d.gHi);
        [X_V, y_V] = stack2(d.vLo, d.vHi);

        % Cross-set: train on one set, test on the other (shuffle TRAIN labels only for the null)
        obs_c  = crossgroup_acc(X_G, y_G, X_V, y_V, 2, y_G, y_V);
        null_c = nan(n_perm_sn,1);
        for pi_ = 1:n_perm_sn
            null_c(pi_) = crossgroup_acc(X_G, y_G, X_V, y_V, 2, ...
                y_G(randperm(numel(y_G))), y_V(randperm(numel(y_V))));
        end
        sn_cross_acc(neuron_i,k) = obs_c;
        sn_cross_p(neuron_i,k)   = (1 + sum(null_c >= obs_c)) / (1 + n_perm_sn);

        % Within-set (pooled, position still aliased with sequence)
        X_all = [X_G; X_V];
        y_all = [y_G; y_V];
        obs_w  = within_cv_acc(X_all, y_all, 2, k_folds_sn, cv_repeats_sn);
        null_w = nan(n_perm_sn,1);
        for pi_ = 1:n_perm_sn
            null_w(pi_) = within_cv_acc(X_all, y_all(randperm(numel(y_all))), 2, k_folds_sn, cv_repeats_sn);
        end
        sn_within_acc(neuron_i,k) = obs_w;
        sn_within_p(neuron_i,k)   = (1 + sum(null_w >= obs_w)) / (1 + n_perm_sn);
    end
end

%% ---- Step 3: summarise single neurons by area ----

n_used  = nan(n_ctx,2);
k_cross = nan(n_ctx,2);
k_with  = nan(n_ctx,2);
acc_cross_by_area  = cell(n_ctx,2);
acc_within_by_area = cell(n_ctx,2);
p_bin_cross_all    = nan(n_ctx,2);
p_fisher_cross_all = nan(n_ctx,1);

for k = 1:n_ctx
    fprintf('\n=== Single-neuron position decoding: %s after %s, position %d vs %d (chance 0.5)%s ===\n', ...
        contexts(k).letter, contexts(k).prev, contexts(k).pos_lo, contexts(k).pos_hi, ...
        ternary_str(k==1,'  [PRIMARY]','  [secondary]'));
    for a = 1:2
        idx = area_idx{a};
        ok  = idx(~isnan(sn_cross_acc(idx,k)));
        n_used(k,a)  = numel(ok);
        k_cross(k,a) = sum(sn_cross_p(ok,k)  < alpha_sn);
        k_with(k,a)  = sum(sn_within_p(ok,k) < alpha_sn);
        acc_cross_by_area{k,a}  = sn_cross_acc(ok,k);
        acc_within_by_area{k,a} = sn_within_acc(ok,k);

        p_bin_c = 1 - binocdf(k_cross(k,a) - 1, n_used(k,a), alpha_sn);
        p_bin_w = 1 - binocdf(k_with(k,a)  - 1, n_used(k,a), alpha_sn);
        p_bin_cross_all(k,a) = p_bin_c;

        cells_med = median(squeeze(sn_cell_n(ok,:,k)), 1, 'omitnan');
        fprintf('\n%s: %d / %d neurons usable (>= %d trials in each of the 4 cells); median trials per cell [gLo gHi vLo vHi] = %s\n', ...
            area_names{a}, n_used(k,a), numel(idx), min_trials_sn, mat2str(cells_med));
        fprintf('  CROSS-set  : median acc = %.3f | %d neurons p<%.2f (%.1f%%, base rate %.0f%%), binomial p = %.4f\n', ...
            median(acc_cross_by_area{k,a}), k_cross(k,a), alpha_sn, 100*k_cross(k,a)/n_used(k,a), 100*alpha_sn, p_bin_c);
        fprintf('  WITHIN-set : median acc = %.3f | %d neurons p<%.2f (%.1f%%, base rate %.0f%%), binomial p = %.4f\n', ...
            median(acc_within_by_area{k,a}), k_with(k,a), alpha_sn, 100*k_with(k,a)/n_used(k,a), 100*alpha_sn, p_bin_w);
    end

    [~, p_fc] = fishertest([k_cross(k,2) n_used(k,2)-k_cross(k,2); k_cross(k,1) n_used(k,1)-k_cross(k,1)]);
    [~, p_fw] = fishertest([k_with(k,2)  n_used(k,2)-k_with(k,2);  k_with(k,1)  n_used(k,1)-k_with(k,1)]);
    p_fisher_cross_all(k) = p_fc;
    fprintf('\n--- Frontal vs auditory, single neurons ---\n');
    fprintf('CROSS-set  fraction significant: frontal %.1f%% vs auditory %.1f%%, Fisher exact p = %.4f\n', ...
        100*k_cross(k,2)/n_used(k,2), 100*k_cross(k,1)/n_used(k,1), p_fc);
    fprintf('WITHIN-set fraction significant: frontal %.1f%% vs auditory %.1f%%, Fisher exact p = %.4f\n', ...
        100*k_with(k,2)/n_used(k,2), 100*k_with(k,1)/n_used(k,1), p_fw);

    aF = acc_cross_by_area{k,2}; aA = acc_cross_by_area{k,1};
    if numel(aF) > 1 && numel(aA) > 1
        fprintf('CROSS-set per-neuron accuracy, frontal (median %.3f) vs auditory (median %.3f): rank-sum p = %.4f\n', ...
            median(aF), median(aA), ranksum(aF, aA));
        mF = nan(n_neuron_boot,1); mA = nan(n_neuron_boot,1);
        for b = 1:n_neuron_boot
            mF(b) = median(aF(randi(numel(aF), numel(aF), 1)));
            mA(b) = median(aA(randi(numel(aA), numel(aA), 1)));
        end
        fprintf('Neuron-resampling bootstrap median difference (frontal - auditory):\n');
        bootstrap_compare(mF, mA);
    end
end

%% ---- Step 4: PCA-free pseudo-population decoding ----

cls = {'lo','hi'};
pop = struct();
for k = 1:n_ctx

    trG = repmat(struct('lo',[],'hi',[]), n_neurons, 1);
    trV = repmat(struct('lo',[],'hi',[]), n_neurons, 1);
    trHi = repmat(struct('G',[],'V',[]), n_neurons, 1);   % context-shift controls
    trLo = repmat(struct('G',[],'V',[]), n_neurons, 1);
    for neuron_i = 1:n_neurons
        d = cdat{neuron_i,k};
        if isempty(d); continue; end
        trG(neuron_i).lo = d.gLo;  trG(neuron_i).hi = d.gHi;
        trV(neuron_i).lo = d.vLo;  trV(neuron_i).hi = d.vHi;
        trHi(neuron_i).G = d.gHi;  trHi(neuron_i).V = d.vHi;
        trLo(neuron_i).G = d.gLo;  trLo(neuron_i).V = d.vLo;
    end

    for a = 1:2
        fprintf('\n=== Population (PCA-free), %s | context %d: %s after %s, pos %d vs %d%s ===\n', ...
            area_names{a}, k, contexts(k).letter, contexts(k).prev, contexts(k).pos_lo, contexts(k).pos_hi, ...
            ternary_str(k==1,' [PRIMARY]',' [secondary]'));

        [b_c, p_c, n_c] = crossgroup_pop_decode(trG, trV, cls, area_idx{a}, min_trials_pop, ...
            n_pseudo_trial, n_time_bins, n_boot, n_perm);
        [b_g, p_g, n_g] = splithalf_pop_decode(trG, cls, area_idx{a}, 2*min_trials_pop, ...
            n_pseudo_trial, n_time_bins, n_boot, n_perm);
        [b_v, p_v, n_v] = splithalf_pop_decode(trV, cls, area_idx{a}, 2*min_trials_pop, ...
            n_pseudo_trial, n_time_bins, n_boot, n_perm);
        [b_h, p_h, n_h] = splithalf_pop_decode(trHi, {'G','V'}, area_idx{a}, 2*min_trials_pop, ...
            n_pseudo_trial, n_time_bins, n_boot, n_perm);
        [b_l, p_l, n_l] = splithalf_pop_decode(trLo, {'G','V'}, area_idx{a}, 2*min_trials_pop, ...
            n_pseudo_trial, n_time_bins, n_boot, n_perm);

        pop(k,a).cross_boot = b_c; pop(k,a).cross_perm = p_c; pop(k,a).n_cross = n_c;

        report_pop('CROSS-set position  (train gram-equiv -> test post-viol, and reverse)', b_c, p_c, n_c, numel(area_idx{a}));
        report_pop('WITHIN gram-equiv set position (aliased with sequence)               ', b_g, p_g, n_g, numel(area_idx{a}));
        report_pop('WITHIN post-violation set position (aliased with sequence)           ', b_v, p_v, n_v, numel(area_idx{a}));
        report_pop(sprintf('CONTEXT SHIFT at position %d (gram-equiv vs post-viol, same letter|pred)', contexts(k).pos_hi), b_h, p_h, n_h, numel(area_idx{a}));
        report_pop(sprintf('CONTEXT SHIFT at position %d (gram-equiv vs post-viol, same letter|pred)', contexts(k).pos_lo), b_l, p_l, n_l, numel(area_idx{a}));
    end

    fprintf('\n--- Population CROSS-set position decoding, context %d: frontal vs auditory ---\n', k);
    bootstrap_compare(pop(k,2).cross_boot, pop(k,1).cross_boot);
end

%% ---- Step 5: evaluate the pre-specified rule on the PRIMARY context ----

k = 1;
accF  = median(pop(k,2).cross_boot, 'omitnan');
accA  = median(pop(k,1).cross_boot, 'omitnan');
permF = prctile(pop(k,2).cross_perm, 95);
dd    = pop(k,2).cross_boot - pop(k,1).cross_boot;
nd    = numel(dd);
p_diff = max(2*min(mean(dd <= 0), mean(dd >= 0)), 1/nd);

crit1 = accF > permF;
crit2 = (accF > accA) && (p_diff < 0.05);
fracF = k_cross(k,2)/n_used(k,2);
fracA = k_cross(k,1)/n_used(k,1);
crit3 = (p_bin_cross_all(k,2) < 0.05) && (fracF > fracA) && (p_fisher_cross_all(k) < 0.05);

fprintf('\n================ PRE-SPECIFIED RULE, PRIMARY CONTEXT (%s after %s, pos %d vs %d) ================\n', ...
    contexts(k).letter, contexts(k).prev, contexts(k).pos_lo, contexts(k).pos_hi);
fprintf('(i)   frontal population cross-set accuracy %.3f vs permutation 95th pct %.3f            : %s\n', accF, permF, pass_str(crit1));
fprintf('(ii)  frontal %.3f vs auditory %.3f, bootstrap two-sided p = %.4f                          : %s\n', accF, accA, p_diff, pass_str(crit2));
fprintf('(iii) frontal significant-neuron fraction %.1f%% (binomial p = %.4f) vs auditory %.1f%% (Fisher p = %.4f) : %s\n', ...
    100*fracF, p_bin_cross_all(k,2), 100*fracA, p_fisher_cross_all(k), pass_str(crit3));
if crit1 && crit2 && crit3
    fprintf('==> SUPPORTED: a position/elapsed-time signal that generalises across non-overlapping sequences is stronger in frontal than auditory.\n');
else
    fprintf('==> NOT SUPPORTED under the rule fixed in advance.\n');
end

%% ---- Step 6: plots (standard MATLAB plotting) ----

cols = [0.2 0.4 0.7; 0.8 0.3 0.2];   % auditory, frontal

figure('Renderer','painters','Position',[100 100 900 400]);
subplot(1,2,1); hold on
for k = 1:n_ctx
    for a = 1:2
        x = (k-1)*3 + a;
        m = median(pop(k,a).cross_boot,'omitnan');
        bar(x, m, 0.8, 'FaceColor', cols(a,:));
        errorbar(x, m, m - prctile(pop(k,a).cross_boot,2.5), prctile(pop(k,a).cross_boot,97.5) - m, 'k', 'LineWidth', 1.2);
        plot(x + [-0.4 0.4], prctile(pop(k,a).cross_perm,95)*[1 1], 'k--');
    end
end
plot([0.3 (n_ctx-1)*3+2.7], [0.5 0.5], 'r--');
set(gca,'XTick',[1.5 4.5],'XTickLabel',{'C|F  pos4 v 5 (primary)','G|C  pos4 v 5'}); ylabel('Cross-set accuracy'); ylim([0 1]); box off
title('Population, cross-set (dashed black = perm 95th)')

subplot(1,2,2); hold on
for k = 1:n_ctx
    for a = 1:2
        x = (k-1)*3 + a;
        bar(x, 100*k_cross(k,a)/n_used(k,a), 0.8, 'FaceColor', cols(a,:));
    end
end
plot([0.3 (n_ctx-1)*3+2.7], [100*alpha_sn 100*alpha_sn], 'r--');
set(gca,'XTick',[1.5 4.5],'XTickLabel',{'C|F  pos4 v 5 (primary)','G|C  pos4 v 5'}); ylabel('% neurons significant (cross-set)'); box off
title('Single neurons (red = 5% base rate); blue = auditory, red-orange = frontal')


%% ================= Local functions =================

function [pre_conds, post_conds] = context_cells(cond_letters, first_viol_pos, letter, prev, pos)
% Conditions containing `letter` preceded by `prev` at position `pos`, split into
% those whose history up to and including that element is grammatical (pre-
% violation or grammatical sequence) and those that follow a violation.
pre_conds = []; post_conds = [];
for c = 1:size(cond_letters,1)
    if strcmp(cond_letters{c,pos}, letter) && strcmp(cond_letters{c,pos-1}, prev)
        if c <= 8 || first_viol_pos(c) > pos
            pre_conds = [pre_conds, c]; %#ok<AGROW>
        else
            post_conds = [post_conds, c]; %#ok<AGROW>
        end
    end
end
end


function s = ternary_str(cond, a, b)
if cond; s = a; else; s = b; end
end


function s = pass_str(tf)
if tf; s = 'PASS'; else; s = 'FAIL'; end
end


function report_pop(label, boot_acc, perm_acc, n_used, n_total)
obs = median(boot_acc,'omitnan');
p   = (1 + sum(perm_acc >= obs)) / (1 + sum(~isnan(perm_acc)));
fprintf('%s\n   neurons usable %d / %d | acc median = %.3f, 95%% CI [%.3f, %.3f] | perm null median %.3f, 95th pct %.3f | perm p = %.4f\n', ...
    label, n_used, n_total, obs, prctile(boot_acc,2.5), prctile(boot_acc,97.5), ...
    median(perm_acc,'omitnan'), prctile(perm_acc,95), p);
end


function [X, y] = stack2(lo, hi)
X = [lo; hi];
y = [ones(size(lo,1),1); 2*ones(size(hi,1),1)];
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


function pred = diag_nc_predict(Xtr, ytr, Xte, K)
% Diagonal-covariance nearest-centroid classifier (diagonal LDA, equal priors).
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
% Train on A -> test on B, and vice versa; mean balanced accuracy. For the
% permutation null, pass shuffled yA_train / yB_train and the TRUE y_A / y_B.
pAB = diag_nc_predict(X_A, yA_train, X_B, K);
pBA = diag_nc_predict(X_B, yB_train, X_A, K);
acc = mean([balanced_acc(pAB, y_B, K), balanced_acc(pBA, y_A, K)]);
end


function acc = within_cv_acc(X, y, K, k_folds, repeats)
% Repeated stratified k-fold CV over pooled trials (sequence NOT controlled).
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
% PCA-free pseudo-population cross-set decode: train on set A, test on set B, and vice versa.

n_classes = length(classes);
K = n_classes;

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


function [boot_acc, perm_acc, n_used] = splithalf_pop_decode(trials, classes, neuron_idx, min_trials, ...
    n_pseudo_trial, n_time_bins, n_boot, n_perm)
% PCA-free pseudo-population decode with a leak-free trial split: each
% neuron's trials in each class are split into disjoint train / test halves on
% every iteration, then pseudo-trials are built from each half separately.

K = numel(classes);
usable = true(numel(neuron_idx),1);
pool   = cell(numel(neuron_idx), K);
for ni = 1:numel(neuron_idx)
    for ci = 1:K
        v = trials(neuron_idx(ni)).(classes{ci});
        pool{ni,ci} = v;
        if isempty(v) || size(v,1) < min_trials
            usable(ni) = false;
        end
    end
end
use_idx = find(usable);
n_used  = numel(use_idx);
pool    = pool(use_idx,:);

boot_acc = nan(n_boot,1);
perm_acc = nan(n_perm,1);
if n_used < 2
    return
end

for it = 1:n_boot
    [ptr, pte] = split_pool(pool);
    [Xtr, ytr] = build_pseudo_population(ptr, n_pseudo_trial, n_time_bins);
    [Xte, yte] = build_pseudo_population(pte, n_pseudo_trial, n_time_bins);
    boot_acc(it) = balanced_acc(diag_nc_predict(Xtr, ytr, Xte, K), yte, K);
end
for it = 1:n_perm
    [ptr, pte] = split_pool(pool);
    [Xtr, ytr] = build_pseudo_population(ptr, n_pseudo_trial, n_time_bins);
    [Xte, yte] = build_pseudo_population(pte, n_pseudo_trial, n_time_bins);
    perm_acc(it) = balanced_acc(diag_nc_predict(Xtr, ytr(randperm(numel(ytr))), Xte, K), yte, K);
end
end


function [ptr, pte] = split_pool(pool)
ptr = cell(size(pool));
pte = cell(size(pool));
for i = 1:size(pool,1)
    for j = 1:size(pool,2)
        v   = pool{i,j};
        n   = size(v,1);
        idx = randperm(n);
        h   = floor(n/2);
        ptr{i,j} = v(idx(1:h),:);
        pte{i,j} = v(idx(h+1:end),:);
    end
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
