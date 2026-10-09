%% Cross-sequence-group generalization test of the PC1-vs-time "linearity" result
%
% position_crossgroup_generalization.m found no generalizable ordinal-
% position decoding above chance in either area once sequence-identity
% could no longer be used as a shortcut. That result only covers the
% manuscript's CATEGORICAL decoding claim (Fig 2C-D). The paper's other,
% structurally different line of evidence for frontal sequential coding
% is the CONTINUOUS trajectory/linearity result (Fig 2A/B, pca_seq_main.m):
% population PC1, built from a per-neuron trace averaged across ALL
% nonviolation trials pooled over all 4 sequences, regressed against
% elapsed time-in-sequence (R^2 = ~0.76 reported for frontal PC1).
%
% Because that trace pools all 4 sequences together per neuron before
% PCA, it isn't confounded by sequence-identity in the same direct way
% classification was -- but it inherits a different risk: if the
% "progressive drift" that PC1 captures is actually driven by mean-level
% differences between specific sequences' acoustic content (which differ
% after position 1) rather than a genuine, sequence-independent temporal/
% ordinal code, then the PC1 axis found by pooling all 4 sequences may
% not be a real population property at all -- it could just be an
% artifact of which particular 4 sequences happen to be in the pool.
%
% This script applies the same structural control used for position
% decoding: split the 4 sequences into the same two disjoint pairs,
%   Set A = sequences {1, 2} (cond_value 1,5,2,6)
%   Set B = sequences {3, 4} (cond_value 3,7,4,8)
% fit PCA (find PC1's axis across neurons) on Set A's pooled trace only,
% then PROJECT Set B's pooled trace onto that same axis (and vice versa).
% If the PC1-vs-time relationship reflects a genuine population property,
% the axis learned from one pair of sequences should still explain
% temporal variance in a completely different pair. If it doesn't
% generalize, the original whole-pool R^2 was likely inflated by
% sequence-specific content rather than reflecting real, generalizable
% temporal/ordinal structure.
%
% Feature extraction, baseline normalization (smooth(...,100), baseline
% -200:0 ms) and the PCA window (pca_window = -100:5:2665, matching
% perform_pca_and_plot.m) are all identical to pca_seq_main.m; the only
% change is that two independent, sequence-disjoint pooled traces are
% built per neuron instead of one. R^2 is invariant to PC sign, so no
% special handling of PCA sign ambiguity is needed.
%
% Assumes the following are already in the workspace: spike_log, dirs,
% auditory_neuron_idx, frontal_neuron_idx.

clear pca_sdf_out_A pca_sdf_out_B n_trials_seq linearity_crossgroup_results

rng(20,'twister')

%% ---- User-set parameters ----
pca_window = -100:5:2665;   % identical to perform_pca_and_plot.m
cols       = 1000 + pca_window;
n_boot     = 300;
nSample    = 500;           % neurons per bootstrap draw, matches pca_seq_main.m
n_pcs      = 3;
min_trials_per_seq = 10;    % matches pca_seq_main.m's own exclusion criterion

seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};
groupset_A_cond = [seq_cond_values{1}, seq_cond_values{2}];   % [1 5 2 6]
groupset_B_cond = [seq_cond_values{3}, seq_cond_values{4}];   % [3 7 4 8]

%% ---- Step 1: per-neuron, per-group-set pooled, baseline z-scored, smoothed traces ----

n_neurons = size(spike_log,1);
pca_sdf_out_A = nan(n_neurons, 6001);
pca_sdf_out_B = nan(n_neurons, 6001);
n_trials_seq  = nan(n_neurons, 4);

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

    all_sdf = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    cond_value_all = event_table.cond_value(nonviol_mask);

    baseline_window = 800:1000;
    baseline_mu_fr  = nanmean(nanmean(all_sdf(:,baseline_window)));
    baseline_std_fr = nanstd(nanmean(all_sdf(:,baseline_window)));

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue
    end

    for seq = 1:4
        n_trials_seq(neuron_i,seq) = sum(ismember(cond_value_all, seq_cond_values{seq}));
    end

    trial_idx_A = ismember(cond_value_all, groupset_A_cond);
    trial_idx_B = ismember(cond_value_all, groupset_B_cond);

    if any(trial_idx_A)
        pca_sdf_out_A(neuron_i,:) = smooth( ...
            (nanmean(all_sdf(trial_idx_A,:)) - baseline_mu_fr) ./ baseline_std_fr, 100);
    end
    if any(trial_idx_B)
        pca_sdf_out_B(neuron_i,:) = smooth( ...
            (nanmean(all_sdf(trial_idx_B,:)) - baseline_mu_fr) ./ baseline_std_fr, 100);
    end
end

% Exclude neurons that fail pca_seq_main.m's own <10-trials-per-sequence criterion
nonvalid = any(n_trials_seq < min_trials_per_seq, 2) | any(isnan(n_trials_seq),2);
pca_sdf_out_A(nonvalid,:) = nan;
pca_sdf_out_B(nonvalid,:) = nan;

fprintf('Neurons passing the >=%d-trials-per-sequence criterion (both group-sets): %d / %d\n', ...
    min_trials_per_seq, sum(~nonvalid), n_neurons);

%% ---- Step 2: bootstrap PCA, same-set and cross-group generalization R^2, per area ----

decoding_problems = {'Auditory', auditory_neuron_idx; 'Frontal', frontal_neuron_idx};

linearity_crossgroup_results = struct();

for prob_i = 1:size(decoding_problems,1)

    label      = decoding_problems{prob_i,1};
    neuron_idx = decoding_problems{prob_i,2};

    fprintf('\n=== %s (PC1-vs-time linearity, cross-sequence-group generalization) ===\n', label);

    r2_sameset = nan(n_boot, n_pcs);
    r2_cross   = nan(n_boot, n_pcs);
    r2_null    = nan(n_boot, n_pcs);

    for boot_i = 1:n_boot

        boot_idx = randsample(neuron_idx, nSample, true);

        A_win = pca_sdf_out_A(boot_idx, cols);
        B_win = pca_sdf_out_B(boot_idx, cols);
        valid = all(~isnan(A_win),2) & all(~isnan(B_win),2);
        A_win = A_win(valid,:);
        B_win = B_win(valid,:);

        if size(A_win,1) < 5
            continue
        end

        Xt_A = A_win';   % time x neurons
        Xt_B = B_win';

        mu_A = mean(Xt_A,1);
        [coeff_A, scores_A] = pca(Xt_A);
        n_keep = min(n_pcs, size(coeff_A,2));

        mu_B = mean(Xt_B,1);
        [coeff_B, scores_B] = pca(Xt_B);

        scores_B_onA = (Xt_B - mu_A) * coeff_A(:,1:n_keep);
        scores_A_onB = (Xt_A - mu_B) * coeff_B(:,1:n_keep);

        % Time-shuffled null: same fitted axes, test-set time order destroyed
        shuf_B = Xt_B(randperm(size(Xt_B,1)),:);
        shuf_A = Xt_A(randperm(size(Xt_A,1)),:);
        scores_shufB_onA = (shuf_B - mu_A) * coeff_A(:,1:n_keep);
        scores_shufA_onB = (shuf_A - mu_B) * coeff_B(:,1:n_keep);

        for pc_i = 1:n_keep
            r2_sameset(boot_i,pc_i) = mean([linear_r2(scores_A(:,pc_i), pca_window(:)), ...
                                             linear_r2(scores_B(:,pc_i), pca_window(:))]);
            r2_cross(boot_i,pc_i)   = mean([linear_r2(scores_B_onA(:,pc_i), pca_window(:)), ...
                                             linear_r2(scores_A_onB(:,pc_i), pca_window(:))]);
            r2_null(boot_i,pc_i)    = mean([linear_r2(scores_shufB_onA(:,pc_i), pca_window(:)), ...
                                             linear_r2(scores_shufA_onB(:,pc_i), pca_window(:))]);
        end
    end

    linearity_crossgroup_results.(label).r2_sameset = r2_sameset;
    linearity_crossgroup_results.(label).r2_cross   = r2_cross;
    linearity_crossgroup_results.(label).r2_null    = r2_null;

    for pc_i = 1:n_pcs
        fprintf('PC%d -- same-set R^2: median = %.3f | cross-group R^2: median = %.3f, 95%% CI [%.3f, %.3f] | shuffled-null: median = %.3f\n', ...
            pc_i, median(r2_sameset(:,pc_i),'omitnan'), median(r2_cross(:,pc_i),'omitnan'), ...
            prctile(r2_cross(:,pc_i),2.5), prctile(r2_cross(:,pc_i),97.5), median(r2_null(:,pc_i),'omitnan'));
    end

    fprintf('PC1 cross-group R^2 vs shuffled-null:\n');
    bootstrap_compare(r2_cross(:,1), r2_null(:,1));

end

%% ---- Step 3: plot (PC1 only) ----

figure('Renderer','painters','Position',[100 100 500 450]); hold on

areas = {'Auditory','Frontal'};
bar_x = 1:4;   % Aud same-set, Aud cross-group, Frontal same-set, Frontal cross-group
bar_median = nan(1,4); err_lo = nan(1,4); err_hi = nan(1,4);

col_i = 1;
for a = 1:2
    same = linearity_crossgroup_results.(areas{a}).r2_sameset(:,1);
    cross = linearity_crossgroup_results.(areas{a}).r2_cross(:,1);
    for series = {same, cross}
        v = series{1};
        bar_median(col_i) = median(v,'omitnan');
        err_lo(col_i) = bar_median(col_i) - prctile(v,2.5);
        err_hi(col_i) = prctile(v,97.5) - bar_median(col_i);
        col_i = col_i + 1;
    end
end

bar(bar_x, bar_median, 0.6, 'FaceColor',[0.4 0.6 0.8]);
errorbar(bar_x, bar_median, err_lo, err_hi, 'k', 'LineStyle','none','LineWidth',1.2)
set(gca,'XTick',bar_x,'XTickLabel',{'Aud same-set','Aud cross-group','Frontal same-set','Frontal cross-group'},'XTickLabelRotation',20)
ylabel('PC1 R^2 (vs. elapsed time-in-sequence)')
ylim([0 1]); box off
title('PC1-vs-time linearity: within one sequence-pair vs. generalizing across pairs')

%% ---- Step 4: frontal vs auditory, cross-group PC1 R^2 ----

fprintf('\nCross-group PC1 R^2: frontal vs auditory\n');
bootstrap_compare(linearity_crossgroup_results.Frontal.r2_cross(:,1), linearity_crossgroup_results.Auditory.r2_cross(:,1));


%% ================= Local functions =================

function r2 = linear_r2(y, x)
% Ordinary R^2 of a simple linear fit of y on x -- matches
% fitlm(x,y).Rsquared.Ordinary used in pca_seq_main.m. R^2 is invariant
% to the sign of y, so PCA sign ambiguity across independent fits needs
% no special handling here.
y = y(:); x = x(:);
p = polyfit(x, y, 1);
yhat = polyval(p, x);
ss_res = sum((y - yhat).^2);
ss_tot = sum((y - mean(y)).^2);
r2 = 1 - ss_res/ss_tot;
end
