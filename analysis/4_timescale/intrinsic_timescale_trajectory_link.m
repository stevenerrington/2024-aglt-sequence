%% Intrinsic timescale linked to single-neuron trajectory contribution and response direction
%
% Motivation: a reviewer flagged that the intrinsic-timescale analysis
% (Fig 4) sits disconnected from the rest of the paper -- as presented,
% "frontal has a longer intrinsic timescale than auditory cortex" mostly
% reproduces the general sensory-to-frontal timescale hierarchy already
% established in the literature (Murray et al. 2014; Cirelli/Wasmuht/
% Zeisler-style analyses), and doesn't obviously connect to this paper's
% own population-trajectory result. This script tests two single-neuron
% properties that COULD link the two:
%
%   (A) Trajectory contribution: each neuron's |loading| on frontal/
%       auditory PC1 (from a non-bootstrapped, full-population PCA,
%       matching pca_seq_main.m's method) as a continuous measure of how
%       strongly that neuron's SDF shape covaries with the population
%       trajectory. Tests whether neurons that contribute more to the
%       trajectory have systematically different intrinsic timescales
%       (correlation, plus a high- vs low-contribution median split).
%
%   (B) Response direction: facilitate vs suppress, from the Fig 3 GLM/
%       MDS clustering (glm_element_clustering.m, the live script).
%       facilitate = clusters [3 9 1 6 13 7 5] (facilitated_clusters +
%       facilitated_ramping), suppress = clusters [10 11 8 4 12 2]
%       (suppressed_clusters + suppressed_ramping) -- i.e. the existing
%       facilitate/suppress grouping already used for Fig 3, pooling
%       across ramping/non-ramping shape (that axis is tested separately
%       in intrinsic_timescale_ramp_vs_nonramp.m).
%
% IMPORTANT FRAMING NOTE (see accompanying message): given the
% pca_ramp_structure_test.m finding that frontal's PC1 ramp is smooth and
% continuous rather than element-locked, "contribution to the PC1
% trajectory" here should NOT be described as "conveys ordinal position
% information" -- that specific claim did not survive earlier testing.
% This script tests whether tau relates to trajectory contribution and to
% response direction, which is a more defensible framing.
%
% Requires, already in the workspace:
%   (a) neuron_tau, neuron_r2, neuron_amp, neuron_total_spikes
%       from intrinsic_timescale_per_neuron.m
%   (b) pca_sdf_out, auditory_neuron_idx, frontal_neuron_idx, and
%       perform_pca_and_plot.m on the path, from (at least the first
%       neuron-loop section of) pca_seq_main.m
%   (c) sig_neurons and mds_results.cluster_idx, from having run (at
%       least through the "Group clusters by response type" section of)
%       glm_element_clustering.m
% If any is missing, this stops with an explicit message rather than
% silently computing something wrong.

%% ---- Check dependencies ----

have_tau = exist('neuron_tau','var') && exist('neuron_r2','var') && ...
           exist('neuron_amp','var') && exist('neuron_total_spikes','var');
have_pca = exist('pca_sdf_out','var') && exist('auditory_neuron_idx','var') && ...
           exist('frontal_neuron_idx','var') && exist('perform_pca_and_plot','file');
have_clusters = exist('sig_neurons','var') && exist('mds_results','var') && isfield(mds_results,'cluster_idx');

if ~have_tau
    error(['neuron_tau / neuron_r2 / neuron_amp / neuron_total_spikes not found. ' ...
           'Run intrinsic_timescale_per_neuron.m first.']);
end
if ~have_pca
    error(['pca_sdf_out / auditory_neuron_idx / frontal_neuron_idx not found, or ' ...
           'perform_pca_and_plot.m is not on the path. Run (at least the neuron-loop ' ...
           'section of) pca_seq_main.m first.']);
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

%% ---- Recompute a full-population (non-bootstrapped) PCA per area ----
% NOTE: by the end of pca_seq_main.m, pc_out_auditory / pc_out_frontal in
% the base workspace hold the LAST bootstrap iteration's PCA (a single
% resampled draw of 500 neurons with replacement), not the full,
% unresampled population. Recomputing fresh here avoids depending on that
% leftover loop state.

pca_window = -100:5:2665;

pc_full_auditory = perform_pca_and_plot(auditory_neuron_idx, pca_sdf_out);
pc_full_frontal  = perform_pca_and_plot(frontal_neuron_idx,  pca_sdf_out);

% Recover which original neurons survived the function's internal
% rmmissing() (neurons with <10 trials/sequence were set to NaN upstream
% and are dropped there) so coeff rows can be mapped back to neuron_i.
valid_aud_idx = auditory_neuron_idx(~any(isnan(pca_sdf_out(auditory_neuron_idx, 1000 + pca_window)), 2));
valid_fro_idx = frontal_neuron_idx(~any(isnan(pca_sdf_out(frontal_neuron_idx,  1000 + pca_window)), 2));

assert(numel(valid_aud_idx) == size(pc_full_auditory.obs.coeff,1), ...
    'Auditory neuron/coeff row mismatch -- PCA neuron filtering has changed upstream.');
assert(numel(valid_fro_idx) == size(pc_full_frontal.obs.coeff,1), ...
    'Frontal neuron/coeff row mismatch -- PCA neuron filtering has changed upstream.');

% Per-neuron |PC1 loading| (PCA sign is arbitrary; magnitude is what
% matters here), indexed by original neuron index.
n_neurons_total = size(spike_log,1);
pc1_loading_abs = nan(n_neurons_total,1);
pc1_loading_abs(valid_aud_idx) = abs(pc_full_auditory.obs.coeff(:,1));
pc1_loading_abs(valid_fro_idx) = abs(pc_full_frontal.obs.coeff(:,1));

%% ========================================================================
%  Section A: intrinsic timescale vs. contribution to the PC1 trajectory
%  ========================================================================

fprintf('\n=== Section A: intrinsic timescale vs. contribution to the PC1 trajectory ===\n');
fprintf('(Contribution = |PC1 loading| from a non-bootstrapped, full-population PCA per area.\n');
fprintf(' Loading sign is arbitrary; magnitude reflects how strongly each neuron''s SDF shape\n');
fprintf(' covaries with the population trajectory -- see framing note in the script header.)\n\n');

aud_sel = intersect(auditory_neuron_idx, find(included & ~isnan(pc1_loading_abs)));
fro_sel = intersect(frontal_neuron_idx,  find(included & ~isnan(pc1_loading_abs)));

area_names = {'Auditory','Frontal'};
area_sels  = {aud_sel, fro_sel};

for area_i = 1:2
    sel  = area_sels{area_i};
    name = area_names{area_i};
    [rho, p_corr] = corr(pc1_loading_abs(sel), neuron_tau(sel), 'type','Spearman');
    fprintf('%s: n=%d, Spearman rho(|PC1 loading|, tau) = %.3f, p = %.4f\n', name, numel(sel), rho, p_corr);
end

fprintf('\n--- Median split: high vs low PC1-loading neurons (bootstrap median tau difference) ---\n');
n_boot = 1000;
rng(4,'twister')

hi_sels = cell(1,2); lo_sels = cell(1,2);

for area_i = 1:2
    sel  = area_sels{area_i};
    name = area_names{area_i};

    med_load = median(pc1_loading_abs(sel));
    hi_sel = sel(pc1_loading_abs(sel) >= med_load);
    lo_sel = sel(pc1_loading_abs(sel) <  med_load);
    hi_sels{area_i} = hi_sel; lo_sels{area_i} = lo_sel;

    fprintf('%s: high-loading n=%d (median tau=%.1fms), low-loading n=%d (median tau=%.1fms)\n', ...
        name, numel(hi_sel), median(neuron_tau(hi_sel)), numel(lo_sel), median(neuron_tau(lo_sel)));

    [p_rs,~,~] = ranksum(neuron_tau(hi_sel), neuron_tau(lo_sel));
    fprintf('  Wilcoxon rank-sum p = %.4f\n', p_rs);

    med_hi_boot = nan(n_boot,1); med_lo_boot = nan(n_boot,1);
    for boot_i = 1:n_boot
        med_hi_boot(boot_i) = median(neuron_tau(randsample(hi_sel, numel(hi_sel), true)));
        med_lo_boot(boot_i) = median(neuron_tau(randsample(lo_sel, numel(lo_sel), true)));
    end
    fprintf('  Bootstrap median difference (High loading - Low loading):\n');
    bootstrap_compare(med_hi_boot, med_lo_boot);
end

%% ========================================================================
%  Section B: intrinsic timescale by response direction (facilitate vs suppress)
%  ========================================================================

fprintf('\n=== Section B: intrinsic timescale by response direction (facilitate vs suppress) ===\n');
fprintf('(Direction from the Fig 3 GLM/MDS clustering: facilitate = clusters [3 9 1 6 13 7 5],\n');
fprintf(' suppress = clusters [10 11 8 4 12 2] -- the existing facilitated/suppressed groups\n');
fprintf(' plus their respective ramping clusters, pooled across ramping/non-ramping shape.)\n\n');

facilitate_clusters_all = [3 9 1 6 13 7 5];
suppress_clusters_all   = [10 11 8 4 12 2];

facilitate_idx = sig_neurons(ismember(mds_results.cluster_idx, facilitate_clusters_all));
suppress_idx   = sig_neurons(ismember(mds_results.cluster_idx, suppress_clusters_all));

area_full_idx = {auditory_neuron_idx, frontal_neuron_idx};
fac_sels = cell(1,2); sup_sels = cell(1,2);

for area_i = 1:2
    area_idx = area_full_idx{area_i};
    name     = area_names{area_i};

    fac_sel = intersect(intersect(facilitate_idx, area_idx), find(included));
    sup_sel = intersect(intersect(suppress_idx,   area_idx), find(included));
    fac_sels{area_i} = fac_sel; sup_sels{area_i} = sup_sel;

    fprintf('%s: facilitate n=%d (median tau=%.1fms), suppress n=%d (median tau=%.1fms)\n', ...
        name, numel(fac_sel), median(neuron_tau(fac_sel)), numel(sup_sel), median(neuron_tau(sup_sel)));

    if numel(fac_sel) < 5 || numel(sup_sel) < 5
        fprintf('  WARNING: fewer than 5 neurons in one group after inclusion filtering -- treat with caution.\n');
    end

    [p_rs,~,~] = ranksum(neuron_tau(fac_sel), neuron_tau(sup_sel));
    fprintf('  Wilcoxon rank-sum p = %.4f\n', p_rs);

    med_fac_boot = nan(n_boot,1); med_sup_boot = nan(n_boot,1);
    for boot_i = 1:n_boot
        med_fac_boot(boot_i) = median(neuron_tau(randsample(fac_sel, numel(fac_sel), true)));
        med_sup_boot(boot_i) = median(neuron_tau(randsample(sup_sel, numel(sup_sel), true)));
    end
    fprintf('  Bootstrap median difference (Facilitate - Suppress):\n');
    bootstrap_compare(med_fac_boot, med_sup_boot);
end

%% ---- Plot 1: tau vs |PC1 loading|, per area ----

figure('Renderer','painters','Position',[100 100 800 380]);

area_colors = [0.2 0.4 0.7; 0.8 0.3 0.2];

for area_i = 1:2
    sel  = area_sels{area_i};
    subplot(1,2,area_i); hold on
    scatter(pc1_loading_abs(sel), neuron_tau(sel), 18, area_colors(area_i,:), 'filled', 'MarkerFaceAlpha', 0.45);
    pfit = polyfit(pc1_loading_abs(sel), neuron_tau(sel), 1);
    xfit = linspace(min(pc1_loading_abs(sel)), max(pc1_loading_abs(sel)), 50);
    plot(xfit, polyval(pfit, xfit), 'k-', 'LineWidth', 1.5);
    xlabel('|PC1 loading|'); ylabel('\tau (ms)'); box off; axis square
    title(area_names{area_i});
end
sgtitle('Per-neuron intrinsic timescale vs. contribution to the PC1 trajectory')

%% ---- Plot 2: facilitate vs suppress, both areas ----

figure('Renderer','painters','Position',[100 100 650 400]); hold on

group_names  = {'Aud. facilitate','Aud. suppress','Fro. facilitate','Fro. suppress'};
group_sel    = {fac_sels{1}, sup_sels{1}, fac_sels{2}, sup_sels{2}};
group_colors = [0.2 0.4 0.7; 0.5 0.6 0.8; 0.8 0.3 0.2; 0.9 0.6 0.5];

for g = 1:4
    v = neuron_tau(group_sel{g});
    jitter_x = g + 0.15*(rand(size(v))-0.5);
    scatter(jitter_x, v, 10, group_colors(g,:), 'filled', 'MarkerFaceAlpha', 0.35);
    plot(g + [-0.25 0.25], median(v)*[1 1], 'k-', 'LineWidth', 2);
end

set(gca,'XTick',1:4,'XTickLabel',group_names,'XTickLabelRotation',15)
ylabel('\tau (ms), per neuron'); box off
title('Per-neuron intrinsic timescale by response direction (black line = median)')
