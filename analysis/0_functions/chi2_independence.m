function [chi2_stat, df, p_value, expected] = chi2_independence(observed)
%CHI2_INDEPENDENCE Performs Chi-square test of independence on a contingency table
%
% INPUT:
%   observed - matrix of observed counts (rows = groups, columns = categories)
%
% OUTPUT:
%   chi2_stat - chi-square statistic
%   df        - degrees of freedom
%   p_value   - p-value
%   expected  - expected counts under independence

    % Validate input
    if any(observed(:) < 0) || ~ismatrix(observed)
        error('Input must be a matrix of non-negative counts.');
    end

    % Row totals and column totals
    row_totals = sum(observed, 2);
    col_totals = sum(observed, 1);
    grand_total = sum(observed, 'all');

    % Expected counts under independence
    expected = (row_totals * col_totals) / grand_total;

    % Chi-square statistic
    chi2_stat = sum((observed - expected).^2 ./ expected, 'all');

    % Degrees of freedom
    df = (size(observed,1)-1) * (size(observed,2)-1);

    % p-value
    p_value = 1 - chi2cdf(chi2_stat, df);
end
