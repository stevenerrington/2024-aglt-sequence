%% Intrinsic timescale linked to single-neuron encoding properties (no population PCA)
%
% Supersedes the PC1-loading part of intrinsic_timescale_trajectory_link.m
% (not yet run) at Steven's request: rather than using each neuron's
% |loading| on a population PCA as the "contributes to ordinal position"
% proxy, this uses a genuinely single-unit measure -- each neuron's own
% GLM encoding of ordinal position, from the SAME sound-modulation GLM
% (glm_sound_modulation.m, run by glm_singleunit_analysis.m) already used
% elsewhere in the pipeline to define "significantly modulated" neurons.
% That GLM fits, per neuron, per 10ms-shifted 100ms window:
%   firing_rate ~ sound     + exp_i   (identity model:  glm_beta rows 1:5,  glm_encoding_flag cols 1:5)
%   firing_rate ~ order_pos + exp_i   (position model:  glm_beta rows 6:10, glm_encoding_flag cols 6:10)
% with FDR correction across time and a >=50ms-significant-run criterion
% (encoding_flag). This script pulls out the position (and, for
% completeness given the RSA identity-coding finding, identity) columns
% per neuron and links them to that neuron's intrinsic timescale.
%
%   (A) Position encoding vs not: binary split on
%       any(glm_encoding_flag(:,6:10),2) -- neurons whose firing rate
%       significantly depends on ordinal position vs those that don't --
%       compared for tau, per area. This is the literal single-unit
%       analogue of "neurons that convey ordinal position vs those that
%       do not."
%   (A2) Continuous version: mean |position beta| in the 0-413ms analysis
%       window as a magnitude measure, correlated with tau, per area.
%   (B) Identity encoding vs not: same as (A) but for the identity
%       (sound A/C/D/F/G) regressor -- supplementary, since the RSA
%       regression found identity coding in auditory but not frontal.
%   (C) Response direction (facilitate vs suppress), from the Fig 3 GLM/
%       MDS clustering -- unchanged from intrinsic_timescale_trajectory_
%       link.m, kept here since it was already single-unit-based.
%
% Requires, already in the workspace:
%   (a) neuron_tau, neuron_r2, neuron_amp, neuron_total_spikes
%       from intrinsic_timescale_per_neuron.m
%   (b) glm_encoding_flag, glm_beta, window_time
%       from glm_singleunit_analysis.m
%   (c) sig_neurons and mds_results.cluster_idx, from having run (at
%       least through the "Group clusters by response type" section of)
%       glm_element_clustering.m
% If any is missing, this stops with an explicit message rather than
% silently computing something wrong.

%% ---- Check dependencies ----

have_tau = exist('neuron_tau','var') && exist('neuron_r2','var') && ...
           exist('neuron_amp','var') && exist('neuron_total_spikes','var');
have_glm = exist('glm_encoding_flag','var') && exist('glm_beta','var') && exist('window_time','var');
have_clusters = exist('sig_neurons','var') && exist('mds_results','var') && isfield(mds_results,'cluster_idx');

if ~have_tau
    error(['neuron_tau / neuron_r2 / neuron_amp / neuron_total_spikes not found. ' ...
           'Run intrinsic_timescale_per_neuron.m first.']);
end
if ~have_glm
    error(['glm_encoding_flag / glm_beta / window_time not found. ' ...
           'Run glm_singleunit_analysis.m first.']);
end
if ~have_clusters
    error(['sig_neurons / mds_results.cluster_idx not found. ' ...
           'Run glm_element_clustering.m first (at least through the ' ...
           '"Group clusters by response type" section).']);
end

%% ---- Re-derive the timescale inclusion criteria (>=250 spikes, R^2>0.5, 0<tau<1000ms, A>0) ----

min_spikes = 250;
r2_cut     = 0.5;
tau_lo     = 0;
tau_hi     = 1000;

included = neuron_total_spikes >= min_spikes & neuron_r2 > r2_cut & ...
           neuron_tau > tau_lo & neuron_tau < tau_hi & neuron_amp > 0;

n_boot = 1000;
rng(5,'twister')

%% ---- Extract per-neuron position / identity GLM encoding (significance + magnitude) ----

n_neurons_total = size(spike_log,1);

glm_timewin_local      = window_time(1,:);
analysis_win_idx_local = find(glm_timewin_local >= 0 & glm_timewin_local <= 413);

% glm_encoding_flag columns: 1:5 = identity (sound A/C/D/F/G), 6:10 = order_pos (position 1-5)
pos_glm_sig      = any(glm_encoding_flag(:,6:10), 2);
identity_glm_sig = any(glm_encoding_flag(:,1:5),  2);

pos_glm_effect      = nan(n_neurons_total,1);
identity_glm_effect = nan(n_neurons_total,1);

for neuron_i = 1:n_neurons_total
    beta_i = glm_beta{neuron_i};
    if isempty(beta_i) || size(beta_i,1) < 10
        continue
    end
    pos_glm_effect(neuron_i)      = mean(abs(beta_i(6:10, analysis_win_idx_local)), 'all', 'omitnan');
    identity_glm_effect(neuron_i) = mean(abs(beta_i(1:5, analysis_win_idx_local)), 'all', 'omitnan');
end

area_names    = {'Auditory','Frontal'};
area_full_idx = {auditory_neuron_idx, frontal_neuron_idx};

%% ========================================================================
%  Section A: intrinsic timescale by position encoding (binary, per-neuron GLM)
%  ========================================================================

fprintf('\n=== Section A: intrinsic timescale by position encoding (per-neuron GLM) ===\n');
fprintf('(Position-encoding = any significant order_pos beta, >=50ms run, FDR-corrected --\n');
fprintf(' the same criterion glm_singleunit_analysis.m uses to flag modulated neurons.)\n\n');

pos_sels = cell(1,2); nonpos_sels = cell(1,2);

for area_i = 1:2
    area_idx = area_full_idx{area_i};
    name     = area_names{area_i};

    pos_sel    = intersect(area_idx, find(included & pos_glm_sig));
    nonpos_sel = intersect(area_idx, find(included & ~pos_glm_sig));
    pos_sels{area_i} = pos_sel; nonpos_sels{area_i} = nonpos_sel;

    compare_tau_groups(neuron_tau, pos_sel, nonpos_sel, ...
        [name ' position-encoding'], [name ' non-position-encoding'], n_boot);
end

%% ---- Section A2: continuous position-effect magnitude vs tau ----

fprintf('\n--- Continuous version: |position beta| vs tau ---\n');

for area_i = 1:2
    sel  = intersect(area_full_idx{area_i}, find(included & ~isnan(pos_glm_effect)));
    name = area_names{area_i};
    [rho, p_corr] = corr(pos_glm_effect(sel), neuron_tau(sel), 'type','Spearman');
    fprintf('%s: n=%d, Spearman rho(|position beta|, tau) = %.3f, p = %.4f\n', name, numel(sel), rho, p_corr);
end

%% ========================================================================
%  Section B: intrinsic timescale by identity encoding (supplementary)
%  ========================================================================

fprintf('\n=== Section B (supplementary): intrinsic timescale by identity encoding ===\n');
fprintf('(Identity-encoding = any significant sound-identity beta, same criterion as above.\n');
fprintf(' Included because the RSA regression found identity coding in auditory but not frontal.)\n\n');

id_sels = cell(1,2); nonid_sels = cell(1,2);

for area_i = 1:2
    area_idx = area_full_idx{area_i};
    name     = area_names{area_i};

    id_sel    = intersect(area_idx, find(included & identity_glm_sig));
    nonid_sel = intersect(area_idx, find(included & ~identity_glm_sig));
    id_sels{area_i} = id_sel; nonid_sels{area_i} = nonid_sel;

    compare_tau_groups(neuron_tau, id_sel, nonid_sel, ...
        [name ' identity-encoding'], [name ' non-identity-encoding'], n_boot);
end

%% ========================================================================
%  Section C: intrinsic timescale by response direction (facilitate vs suppress)
%  ========================================================================

fprintf('\n=== Section C: intrinsic timescale by response direction (facilitate vs suppress) ===\n');
fprintf('(Direction from the Fig 3 GLM/MDS clustering: facilitate = clusters [3 9 1 6 13 7 5],\n');
fprintf(' suppress = clusters [10 11 8 4 12 2] -- the existing facilitated/suppressed groups\n');
fprintf(' plus their respective ramping clusters, pooled across ramping/non-ramping shape.)\n\n');

facilitate_clusters_all = [3 9 1 6 13 7 5];
suppress_clusters_all   = [10 11 8 4 12 2];

facilitate_idx = sig_neurons(ismember(mds_results.cluster_idx, facilitate_clusters_all));
suppress_idx   = sig_neurons(ismember(mds_results.cluster_idx, suppress_clusters_all));

fac_sels = cell(1,2); sup_sels = cell(1,2);

for area_i = 1:2
    area_idx = area_full_idx{area_i};
    name     = area_names{area_i};

    fac_sel = intersect(intersect(facilitate_idx, area_idx), find(included));
    sup_sel = intersect(intersect(suppress_idx,   area_idx), find(included));
    fac_sels{area_i} = fac_sel; sup_sels{area_i} = sup_sel;

    compare_tau_groups(neuron_tau, fac_sel, sup_sel, ...
        [name ' facilitate'], [name ' suppress'], n_boot);
end

%% ---- Plot 1: position-encoding vs not, both areas ----

figure('Renderer','painters','Position',[100 100 650 400]); hold on

group_names  = {'Aud. position+','Aud. position-','Fro. position+','Fro. position-'};
group_sel    = {pos_sels{1}, nonpos_sels{1}, pos_sels{2}, nonpos_sels{2}};
group_colors = [0.2 0.4 0.7; 0.5 0.6 0.8; 0.8 0.3 0.2; 0.9 0.6 0.5];

for g = 1:4
    v = neuron_tau(group_sel{g});
    jitter_x = g + 0.15*(rand(size(v))-0.5);
    scatter(jitter_x, v, 10, group_colors(g,:), 'filled', 'MarkerFaceAlpha', 0.35);
    plot(g + [-0.25 0.25], median(v)*[1 1], 'k-', 'LineWidth', 2);
end

set(gca,'XTick',1:4,'XTickLabel',group_names,'XTickLabelRotation',15)
ylabel('\tau (ms), per neuron'); box off
title('Per-neuron intrinsic timescale by position encoding (black line = median)')

%% ---- Plot 2: tau vs |position beta|, per area ----

figure('Renderer','painters','Position',[100 100 800 380]);

area_colors = [0.2 0.4 0.7; 0.8 0.3 0.2];

for area_i = 1:2
    sel  = intersect(area_full_idx{area_i}, find(included & ~isnan(pos_glm_effect)));
    subplot(1,2,area_i); hold on
    scatter(pos_glm_effect(sel), neuron_tau(sel), 18, area_colors(area_i,:), 'filled', 'MarkerFaceAlpha', 0.45);
    pfit = polyfit(pos_glm_effect(sel), neuron_tau(sel), 1);
    xfit = linspace(min(pos_glm_effect(sel)), max(pos_glm_effect(sel)), 50);
    plot(xfit, polyval(pfit, xfit), 'k-', 'LineWidth', 1.5);
    xlabel('|position beta| (a.u.)'); ylabel('\tau (ms)'); box off; axis square
    title(area_names{area_i});
end
sgtitle('Per-neuron intrinsic timescale vs. position-encoding magnitude (single-unit GLM)')

%% ---- Plot 3: facilitate vs suppress, both areas ----

figure('Renderer','painters','Position',[100 100 650 400]); hold on

group_names2  = {'Aud. facilitate','Aud. suppress','Fro. facilitate','Fro. suppress'};
group_sel2    = {fac_sels{1}, sup_sels{1}, fac_sels{2}, sup_sels{2}};

for g = 1:4
    v = neuron_tau(group_sel2{g});
    jitter_x = g + 0.15*(rand(size(v))-0.5);
    scatter(jitter_x, v, 10, group_colors(g,:), 'filled', 'MarkerFaceAlpha', 0.35);
    plot(g + [-0.25 0.25], median(v)*[1 1], 'k-', 'LineWidth', 2);
end

set(gca,'XTick',1:4,'XTickLabel',group_names2,'XTickLabelRotation',15)
ylabel('\tau (ms), per neuron'); box off
title('Per-neuron intrinsic timescale by response direction (black line = median)')

%% ---- Local function: rank-sum + bootstrap median-difference report ----

function compare_tau_groups(neuron_tau, sel_a, sel_b, label_a, label_b, n_boot)
    fprintf('%s n=%d (median tau=%.1fms), %s n=%d (median tau=%.1fms)\n', ...
        label_a, numel(sel_a), median(neuron_tau(sel_a)), label_b, numel(sel_b), median(neuron_tau(sel_b)));

    if numel(sel_a) < 5 || numel(sel_b) < 5
        fprintf('  WARNING: fewer than 5 neurons in one group after inclusion filtering -- treat with caution.\n');
    end

    [p_rs,~,~] = ranksum(neuron_tau(sel_a), neuron_tau(sel_b));
    fprintf('  Wilcoxon rank-sum p = %.4f\n', p_rs);

    med_a_boot = nan(n_boot,1); med_b_boot = nan(n_boot,1);
    for boot_i = 1:n_boot
        med_a_boot(boot_i) = median(neuron_tau(randsample(sel_a, numel(sel_a), true)));
        med_b_boot(boot_i) = median(neuron_tau(randsample(sel_b, numel(sel_b), true)));
    end
    fprintf('  Bootstrap median difference (%s - %s):\n', label_a, label_b);
    bootstrap_compare(med_a_boot, med_b_boot);
end
