function [p_val, ci, diff_dist] = bootstrap_compare(obs_vals, shuf_vals)
%BOOTSTRAP_COMPARE Compare two bootstrap distributions
%
%   [p_val, ci, diff_dist] = bootstrap_compare(obs_vals, shuf_vals)
%
%   obs_vals   : bootstrap samples of observed metric
%   shuf_vals  : bootstrap samples of shuffled metric
%
% Returns:
%   diff_dist  : bootstrap differences (obs - shuf)
%   ci         : 95% confidence interval of the difference
%   p_val      : two-sided bootstrap p-value
%
% Example:
%   [p, ci] = bootstrap_compare(obs_r2.frontal(:,2), shuf_r2.frontal(:,1));

    % Compute bootstrap difference distribution
    diff_dist = obs_vals - shuf_vals;

    % 95% CI on difference
    ci = prctile(diff_dist, [2.5 50.0 97.5]);

    % Two-sided bootstrap p-value
    p_val = 2 * min( mean(diff_dist >= 0), mean(diff_dist <= 0) );

end