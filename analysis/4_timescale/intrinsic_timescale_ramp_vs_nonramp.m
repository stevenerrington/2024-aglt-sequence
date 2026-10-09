%% Intrinsic timescale: (1) auditory vs frontal, (2) frontal ramping vs frontal non-ramping
%
% Comparison (1) is the same auditory-vs-frontal test already run in
% intrinsic_timescale_per_neuron.m, reproduced here for completeness.
%
% Comparison (2) is new: it cross-references each frontal neuron's
% single-unit functional category from the Fig 3 GLM/MDS-clustering
% pipeline (glm_element_clustering.m -- the live script; glm_clustering.m
% and glm_element_clustering_old.m are superseded/dead code and are NOT
% used here) against its per-neuron intrinsic timescale from
% intrinsic_timescale_per_neuron.m. In that pipeline, 13 k-means clusters
% of the GLM-significant neurons' SDF shapes are grouped into 4 response
% types:
%   facilitated_clusters = [3 9 1 6 13 7]   (transient facilitation)
%   suppressed_clusters  = [10 11 8 4 12]   (transient suppression)
%   facilitated_ramping  = 5                (ramping response)
%   suppressed_ramping   = 2                (ramping response)
% "Ramping" here = clusters 5 and 2 combined (the two clusters the
% existing pipeline itself already labels as ramping motifs); "non-
% ramping" = the other 11 clusters (facilitated_clusters + suppressed_
% clusters). This tests whether the longer intrinsic timescale (and the
% continuous, non-element-locked ramp identified in pca_ramp_structure_
% test.m) in frontal cortex is specifically carried by the subset of
% frontal neurons independently classified as ramping at the single-unit
% level, or is a more general property of frontal neurons regardless of
% response shape.
%
% Requires TWO things already in the workspace:
%   (a) neuron_tau, neuron_r2, neuron_amp, neuron_total_spikes (or
%       neuron_tau_table) from intrinsic_timescale_per_neuron.m
%   (b) sig_neurons and mds_results.cluster_idx from having run (at
%       least through the "Group clusters by response type" section of)
%       glm_element_clustering.m
% If either is missing, this script prints what to run first and stops
% rather than silently computing something wrong.
%
% Assumes the following are already in the workspace: auditory_neuron_idx,
% frontal_neuron_idx, plus (a) and (b) above.

%% ---- Check dependencies ----

have_tau = exist('neuron_tau','var') && exist('neuron_r2','var') && ...
           exist('neuron_amp','var') && exist('neuron_total_spikes','var');
have_clusters = exist('sig_neurons','var') && exist('mds_results','var') && isfield(mds_results,'cluster_idx');

if ~have_tau
    error(['neuron_tau / neuron_r2 / neuron_amp / neuron_total_spikes not found. ' ...
           'Run intrinsic_timescale_per_neuron.m first.']);
end
if ~have_clusters
    error(['sig_neurons / mds_results.cluster_idx not found. ' ...
           'Run glm_element_clustering.m first (at least through the ' ...
           '"Group clusters by response type" section).']);
end

%% ---- Re-derive the inclusion criteria (>=250 spikes, R^2>0.5, 0<tau<1000ms, A>0) ----

min_spikes = 250;
r2_cut     = 0.5;
tau_lo     = 0;
tau_hi     = 1000;

included = neuron_total_spikes >= min_spikes & neuron_r2 > r2_cut & ...
           neuron_tau > tau_lo & neuron_tau < tau_hi & neuron_amp > 0;

%% ---- Define the ramp / non-ramp frontal groups from the Fig 3 clustering ----

ramp_clusters    = [5 2];                           % facilitated_ramping, suppressed_ramping
nonramp_clusters = [3 9 1 6 13 7 10 11 8 4 12];      % facilitated_clusters + suppressed_clusters

frontal_ramp_idx    = intersect(sig_neurons(ismember(mds_results.cluster_idx, ramp_clusters)), frontal_neuron_idx);
frontal_nonramp_idx = intersect(sig_neurons(ismember(mds_results.cluster_idx, nonramp_clusters)), frontal_neuron_idx);

fprintf('Frontal neurons classified as ramping (Fig 3 clusters 5+2): %d\n', numel(frontal_ramp_idx));
fprintf('Frontal neurons classified as non-ramping (Fig 3 clusters 3,9,1,6,13,7,10,11,8,4,12): %d\n', numel(frontal_nonramp_idx));

%% ---- Apply the timescale inclusion criteria on top ----

aud_sel      = intersect(auditory_neuron_idx, find(included));
fro_sel      = intersect(frontal_neuron_idx,  find(included));
fro_ramp_sel = intersect(frontal_ramp_idx,    find(included));
fro_nonramp_sel = intersect(frontal_nonramp_idx, find(included));

fprintf('\nAfter applying timescale inclusion criteria (>=%d spikes, R^2>%.1f, %d<tau<%d, A>0):\n', ...
    min_spikes, r2_cut, tau_lo, tau_hi);
fprintf('  Auditory (all):        n = %d\n', numel(aud_sel));
fprintf('  Frontal (all):         n = %d\n', numel(fro_sel));
fprintf('  Frontal, ramping:      n = %d (%.0f%% of classified ramping neurons)\n', ...
    numel(fro_ramp_sel), 100*numel(fro_ramp_sel)/max(numel(frontal_ramp_idx),1));
fprintf('  Frontal, non-ramping:  n = %d (%.0f%% of classified non-ramping neurons)\n', ...
    numel(fro_nonramp_sel), 100*numel(fro_nonramp_sel)/max(numel(frontal_nonramp_idx),1));

%% ---- Comparison 1: auditory vs frontal ----

fprintf('\n=== (1) Auditory vs Frontal ===\n');
fprintf('Auditory median tau = %.1f ms (n=%d) | Frontal median tau = %.1f ms (n=%d)\n', ...
    median(neuron_tau(aud_sel)), numel(aud_sel), median(neuron_tau(fro_sel)), numel(fro_sel));
[p_ranksum_1,~,~] = ranksum(neuron_tau(fro_sel), neuron_tau(aud_sel));
fprintf('Wilcoxon rank-sum p = %.4f\n', p_ranksum_1);

n_boot = 1000;
rng(2,'twister')
med_aud_boot = nan(n_boot,1);
med_fro_boot = nan(n_boot,1);
for boot_i = 1:n_boot
    med_aud_boot(boot_i) = median(neuron_tau(randsample(aud_sel, numel(aud_sel), true)));
    med_fro_boot(boot_i) = median(neuron_tau(randsample(fro_sel, numel(fro_sel), true)));
end
fprintf('Bootstrap median difference (Frontal - Auditory):\n');
bootstrap_compare(med_fro_boot, med_aud_boot);

%% ---- Comparison 2: frontal non-ramping vs frontal ramping ----

fprintf('\n=== (2) Frontal non-ramping vs Frontal ramping ===\n');
fprintf('Non-ramping median tau = %.1f ms (n=%d) | Ramping median tau = %.1f ms (n=%d)\n', ...
    median(neuron_tau(fro_nonramp_sel)), numel(fro_nonramp_sel), ...
    median(neuron_tau(fro_ramp_sel)), numel(fro_ramp_sel));

if numel(fro_ramp_sel) < 5 || numel(fro_nonramp_sel) < 5
    fprintf('WARNING: fewer than 5 neurons in one group after inclusion filtering -- treat this comparison with caution.\n');
end

[p_ranksum_2,~,~] = ranksum(neuron_tau(fro_ramp_sel), neuron_tau(fro_nonramp_sel));
fprintf('Wilcoxon rank-sum p = %.4f\n', p_ranksum_2);

med_nonramp_boot = nan(n_boot,1);
med_ramp_boot    = nan(n_boot,1);
for boot_i = 1:n_boot
    med_nonramp_boot(boot_i) = median(neuron_tau(randsample(fro_nonramp_sel, numel(fro_nonramp_sel), true)));
    med_ramp_boot(boot_i)    = median(neuron_tau(randsample(fro_ramp_sel, numel(fro_ramp_sel), true)));
end
fprintf('Bootstrap median difference (Ramping - Non-ramping):\n');
bootstrap_compare(med_ramp_boot, med_nonramp_boot);

%% ---- Plot: all four groups ----

figure('Renderer','painters','Position',[100 100 650 400]); hold on

group_names   = {'Auditory','Frontal (all)','Frontal non-ramp','Frontal ramp'};
group_sel     = {aud_sel, fro_sel, fro_nonramp_sel, fro_ramp_sel};
group_colors  = [0.2 0.4 0.7; 0.8 0.3 0.2; 0.9 0.6 0.3; 0.6 0.1 0.1];

for g = 1:4
    v = neuron_tau(group_sel{g});
    jitter_x = g + 0.15*(rand(size(v))-0.5);
    scatter(jitter_x, v, 10, group_colors(g,:), 'filled', 'MarkerFaceAlpha', 0.35);
    plot(g + [-0.25 0.25], median(v)*[1 1], 'k-', 'LineWidth', 2);
end

set(gca,'XTick',1:4,'XTickLabel',group_names,'XTickLabelRotation',15)
ylabel('\tau (ms), per neuron'); box off
title('Per-neuron intrinsic timescale by group (black line = median)')
