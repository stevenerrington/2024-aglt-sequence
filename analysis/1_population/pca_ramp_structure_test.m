%% Is the PC1 "ramp" a continuous, generic time-elapsed signal, or discrete, element-locked steps?
%
% Motivation. pca_seq_linearity_crossgroup.m found a real, cross-sequence-
% group-generalizing relationship between PC1 and elapsed time (frontal
% cross-group R^2 = 0.172, p = 0.013 vs a shuffled-time null; auditory
% R^2 = 0.125, p = 0.053) -- so the trial-averaged population trajectory
% really does ramp with time, robustly. But position_crossgroup_
% generalization.m found chance-level single-trial decoding of ordinal
% position in both areas. These are not actually contradictory: the R^2
% analyses regress a trial-AVERAGED, 100ms-smoothed curve against time --
% essentially noise-free at the level of the mean -- while the decoding
% analyses classify single trials, where trial-to-trial variability can
% easily swamp a real but modest mean-level trend. A clean average ramp
% and near-chance single-trial categorical decoding routinely coexist in
% neural data; they are different statistical targets, not the same claim
% at different rigor levels.
%
% There is a second, more specific possibility worth ruling in or out.
% Because every sequence has the same fixed duration and the same fixed
% time-to-reward, "elapsed time since sequence onset" is, in this design,
% perfectly confounded with "time remaining until reward". A regression
% of PC1 against elapsed time cannot distinguish a signal that represents
% which ordinal element is currently being processed from a generic
% ramp-to-bound / reward-anticipation signal that just tracks the passage
% of time and would be identical regardless of sequence content. The two
% make different predictions about the SHAPE of the trajectory:
%   - genuine, element-locked ordinal coding predicts a staircase: activity
%     should change disproportionately in a short window straddling each
%     element onset (563, 1126, 1688, 2252 ms) and stay comparatively flat
%     in between transitions.
%   - a generic, continuous elapsed-time/anticipation signal predicts a
%     smooth ramp: the local rate of change should look about the same
%     within an element's window as it does at the boundary between
%     elements -- no special structure locked to element onsets at all.
%
% This script tests that directly on the same per-neuron pooled trace
% used for the manuscript's reported R^2 = 0.76 (pca_seq_main.m: baseline
% z-scored, smooth(...,100), averaged across ALL nonviolation trials
% pooled over all 4 sequences), inside the same bootstrap-over-neurons /
% PCA-per-draw procedure as pca_seq_linearity_crossgroup.m. In each
% bootstrap draw it computes: (1) the ordinary linear R^2 of PC1 against
% time (a sanity-check reproduction of the existing result), (2) a
% 6-level step/categorical R^2 (pre-onset + one bin per element) --
% which, having more free parameters, can only be >= the linear R^2 --
% and takes the GAP between them as a descriptive measure of how much
% extra, non-straight-line structure is present, and (3) the local slope
% of PC1 vs time in narrow windows straddling each element boundary,
% compared against the local slope in the central 60% of each element's
% own window (comparison performed with bootstrap_compare across the 300
% matched bootstrap draws).
%
% Assumes the following are already in the workspace: spike_log, dirs,
% auditory_neuron_idx, frontal_neuron_idx.

clear pca_sdf_out n_trials_seq ramp_structure_results

rng(21,'twister')

%% ---- User-set parameters ----
pca_window   = -100:5:2665;   % identical to perform_pca_and_plot.m / pca_seq_main.m
cols         = 1000 + pca_window;
n_boot       = 300;
nSample      = 500;           % neurons per bootstrap draw, matches pca_seq_main.m
min_trials_per_seq = 10;      % matches pca_seq_main.m's own exclusion criterion

seq_cond_values   = {[1 5],[2 6],[3 7],[4 8]};
element_bin_edges = [-100 0 563 1126 1688 2252 2665];   % pre, elem1..elem5
boundary_times    = [563 1126 1688 2252];               % internal element-to-element transitions
boundary_halfwidth = 40;      % +/- ms around each boundary used for "transition" slope
mid_fraction_margin = 0.2;    % exclude the outer 20% of each element window; use the central 60% for "mid-element" slope

%% ---- Step 1: per-neuron trace pooled across ALL 4 sequences (matches pca_seq_main.m exactly) ----

n_neurons = size(spike_log,1);
pca_sdf_out  = nan(n_neurons, 6001);
n_trials_seq = nan(n_neurons, 4);

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

    pca_sdf_out(neuron_i,:) = smooth( ...
        (nanmean(all_sdf) - baseline_mu_fr) ./ baseline_std_fr, 100);
end

nonvalid = any(n_trials_seq < min_trials_per_seq, 2) | any(isnan(n_trials_seq),2);
pca_sdf_out(nonvalid,:) = nan;

fprintf('Neurons passing >=%d-trials-per-sequence criterion (all 4 pooled): %d / %d\n', ...
    min_trials_per_seq, sum(~nonvalid), n_neurons);

%% ---- Step 2: precompute fixed time-axis masks (same for every bootstrap draw) ----

bin_idx = discretize(pca_window, element_bin_edges);   % 1 (pre) .. 6 (elem5), NaN if outside range

boundary_mask = false(size(pca_window));
for b = 1:numel(boundary_times)
    boundary_mask = boundary_mask | (abs(pca_window - boundary_times(b)) <= boundary_halfwidth);
end

elem_edges = element_bin_edges(2:end);   % [0 563 1126 1688 2252 2665]
mid_mask = false(size(pca_window));
for k = 1:5
    a = elem_edges(k); b = elem_edges(k+1);
    lo = a + mid_fraction_margin*(b-a);
    hi = b - mid_fraction_margin*(b-a);
    mid_mask = mid_mask | (pca_window >= lo & pca_window < hi);
end

%% ---- Step 3: bootstrap PCA per area; linear vs step R^2, boundary vs mid-element local slope ----

decoding_problems = {'Auditory', auditory_neuron_idx; 'Frontal', frontal_neuron_idx};

ramp_structure_results = struct();

for prob_i = 1:size(decoding_problems,1)

    label      = decoding_problems{prob_i,1};
    neuron_idx = decoding_problems{prob_i,2};

    fprintf('\n=== %s (continuous-ramp vs. element-locked-step test) ===\n', label);

    r2_linear_all      = nan(n_boot,1);
    r2_step_all        = nan(n_boot,1);
    extra_r2_all       = nan(n_boot,1);
    slope_boundary_all = nan(n_boot,1);
    slope_mid_all      = nan(n_boot,1);
    pc1_traj_all       = nan(n_boot, numel(pca_window));

    for boot_i = 1:n_boot

        boot_idx = randsample(neuron_idx, nSample, true);
        win = pca_sdf_out(boot_idx, cols);
        valid = all(~isnan(win),2);
        win = win(valid,:);

        if size(win,1) < 5
            continue
        end

        Xt = win';                     % time x neurons
        [~, scores] = pca(Xt);
        pc1 = scores(:,1);             % length(pca_window) x 1

        % Fix PC1 sign so slopes are comparable across bootstrap draws
        % (R^2 measures below are sign-invariant, but signed local slopes are not)
        p_overall = polyfit(pca_window(:), pc1, 1);
        if p_overall(1) < 0
            pc1 = -pc1;
        end

        r2_linear_all(boot_i) = linear_r2(pc1, pca_window(:));

        grand_mean = mean(pc1);
        ss_tot = sum((pc1 - grand_mean).^2);
        ss_between = 0;
        for bi = 1:6
            sel = bin_idx == bi;
            if any(sel)
                bin_mean = mean(pc1(sel));
                ss_between = ss_between + sum(sel)*(bin_mean - grand_mean)^2;
            end
        end
        r2_step_all(boot_i)  = ss_between / ss_tot;
        extra_r2_all(boot_i) = r2_step_all(boot_i) - r2_linear_all(boot_i);

        p_boundary = polyfit(pca_window(boundary_mask)', pc1(boundary_mask), 1);
        p_mid      = polyfit(pca_window(mid_mask)',      pc1(mid_mask),      1);
        slope_boundary_all(boot_i) = p_boundary(1);
        slope_mid_all(boot_i)      = p_mid(1);

        pc1_traj_all(boot_i,:) = pc1';
    end

    ramp_structure_results.(label).r2_linear      = r2_linear_all;
    ramp_structure_results.(label).r2_step        = r2_step_all;
    ramp_structure_results.(label).extra_r2       = extra_r2_all;
    ramp_structure_results.(label).slope_boundary = slope_boundary_all;
    ramp_structure_results.(label).slope_mid      = slope_mid_all;
    ramp_structure_results.(label).pc1_traj       = pc1_traj_all;

    fprintf('Linear R^2 (PC1 vs time): median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        median(r2_linear_all,'omitnan'), prctile(r2_linear_all,2.5), prctile(r2_linear_all,97.5));
    fprintf('Step (6-level) R^2: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        median(r2_step_all,'omitnan'), prctile(r2_step_all,2.5), prctile(r2_step_all,97.5));
    fprintf('Extra R^2 from step beyond linear: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        median(extra_r2_all,'omitnan'), prctile(extra_r2_all,2.5), prctile(extra_r2_all,97.5));
    fprintf('Local slope (PC1 units / ms) -- boundary-window median = %.5f | mid-element-window median = %.5f\n', ...
        median(slope_boundary_all,'omitnan'), median(slope_mid_all,'omitnan'));
    fprintf('Boundary-window slope vs mid-element-window slope:\n');
    bootstrap_compare(slope_boundary_all, slope_mid_all);

end

%% ---- Step 4: frontal vs auditory, on the structure diagnostics themselves ----

fprintf('\nExtra R^2 (step beyond linear): frontal vs auditory\n');
bootstrap_compare(ramp_structure_results.Frontal.extra_r2, ramp_structure_results.Auditory.extra_r2);

fprintf('\nBoundary-vs-mid slope ratio: frontal vs auditory\n');
bootstrap_compare(ramp_structure_results.Frontal.slope_boundary - ramp_structure_results.Frontal.slope_mid, ...
                   ramp_structure_results.Auditory.slope_boundary - ramp_structure_results.Auditory.slope_mid);

%% ---- Step 5: plot median PC1 trajectory per area, with element boundaries marked ----

figure('Renderer','painters','Position',[100 100 650 450]); hold on

area_colors = [0.2 0.4 0.7; 0.8 0.3 0.2];   % Auditory, Frontal
areas = {'Auditory','Frontal'};

y_lims = [Inf -Inf];
for a = 1:2
    traj = ramp_structure_results.(areas{a}).pc1_traj;
    med_traj = median(traj,1,'omitnan');
    lo_traj  = prctile(traj,2.5,1);
    hi_traj  = prctile(traj,97.5,1);
    y_lims(1) = min(y_lims(1), min(lo_traj));
    y_lims(2) = max(y_lims(2), max(hi_traj));

    patch([pca_window, fliplr(pca_window)], [lo_traj, fliplr(hi_traj)], area_colors(a,:), ...
        'FaceAlpha', 0.15, 'EdgeColor','none');
    plot(pca_window, med_traj, 'Color', area_colors(a,:), 'LineWidth', 1.8, 'DisplayName', areas{a});
end

for b = 1:numel(boundary_times)
    xline(boundary_times(b), 'k--', 'LineWidth', 1);
end
xline(0, 'k-', 'LineWidth', 1);

xlabel('Time from sequence onset (ms)')
ylabel('PC1 (sign-corrected, a.u.)')
title('PC1 trajectory: smooth ramp or element-locked steps? (dashed = element boundaries)')
legend(areas, 'Location','northwest')
box off
ylim(y_lims)


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
