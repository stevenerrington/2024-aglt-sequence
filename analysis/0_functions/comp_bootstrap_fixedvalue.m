function comp_bootstrap_fixedvalue(bootValues, refValue)
    % reportBootstrap - Computes 95% CI and two-sided bootstrap p-value
    % bootValues : vector of bootstrapped values (can contain NaN)
    % refValue   : value to test against
    %
    % Example output: "CI [512.3 to 579.1 ms], p = 0.024"
    
    % Remove NaNs for percentile calculation
    cleanBoot = bootValues(~isnan(bootValues));
    
    % 95% confidence interval
    CI = prctile(cleanBoot, [2.5 50.0 97.5]);
    
    % Center of bootstrap distribution using nanmean
    center = nanmean(cleanBoot);
    
    % Two-sided p-value
     p = mean(cleanBoot >= refValue);   % proportion greater than refValue
    
    % Display nicely
    fprintf('CI [%.2f, %.2f, %.2f ms], p = %.3f\n', CI(1), CI(2), CI(3), p);


end
