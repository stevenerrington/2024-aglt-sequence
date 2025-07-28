function [chi2stat, p, stats] = chi2cont(table)
    % Performs chi-squared test of independence on contingency table

    [rows, cols] = size(table);
    total = sum(table(:));
    row_sums = sum(table, 2);
    col_sums = sum(table, 1);
    expected = row_sums * col_sums / total;

    chi2stat = sum((table - expected).^2 ./ expected, 'all');
    df = (rows - 1) * (cols - 1);
    p = 1 - chi2cdf(chi2stat, df);
    
    stats.expected = expected;
    stats.df = df;
    stats.chi2stat = chi2stat;
end
