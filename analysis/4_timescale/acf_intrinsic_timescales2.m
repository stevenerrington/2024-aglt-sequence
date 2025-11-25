% ========================================================================
% Script: bootstrap_acf_timescale.m
%
% Purpose:
%   Bootstrap spike-count autocorrelation functions (ACFs) from baseline
%   spike activity, estimate neuronal timescales (tau) by fitting an
%   exponential decay model, and compare timescales between auditory and
%   frontal neuron populations.
%
% High-Level Workflow:
%   1. Load / cache results:
%        - If cached ACFs (acf_boot_out.mat) do not exist, compute them.
%   2. Per neuron, per bootstrap iteration:
%        - Sample 25 valid trials (with replacement) meeting condition:
%          cond_label == 'nonviol' & rewardOnset_ms finite.
%        - Bin spikes in baseline window [-1000, 0] ms using 50 ms bins.
%        - Compute across-trial Pearson correlations for all bin pairs.
%        - Collapse correlation matrix by absolute lag -> 1D ACF vs lag (ms).
%        - Store ACF in acf_boot_out(neuron, lag, boot).
%   3. Save bootstrapped ACFs to disk (or load if already computed).
%   4. For each area (auditory, frontal):
%        - Select neurons_in = intersection of area idx and nonzero_neurons.
%        - For each bootstrap: average ACF across selected neurons (fit range
%          50–500 ms) and fit exponential decay A*exp(-lag/tau)+C via lsqcurvefit.
%        - Store tau in tau_test(:, area_i).
%   5. Visualize bootstrap tau distributions (histograms).
%   6. Compute frontal–auditory tau difference distribution, 95% CI, and a
%      simple two-tailed p-estimate based on sign of differences.
%   7. Plot tau difference using gramm (jitter + percentile summary).
%
% Dependencies:
%   MATLAB Optimization Toolbox (lsqcurvefit), gramm (for plotting), figuren()
%   utility (custom, assumed in path), Statistics Toolbox (randsample, histcounts).
%
% Author: Steven Errington   Date: 2025-07-21
% ========================================================================

