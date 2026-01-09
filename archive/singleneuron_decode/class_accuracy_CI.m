function [ci_lower, ci_upper] = class_accuracy_CI(X, Y, nFolds, classLabel, nPerm)
% Compute 95% confidence interval for decoding a single class
%
% Inputs:
%   X         - feature matrix (nTrials x nFeatures)
%   Y         - labels (categorical)
%   nFolds    - cross-validation folds
%   classLabel- the class of interest (categorical or string)
%   nPerm     - number of permutations
%
% Outputs:
%   ci_lower, ci_upper - 95% confidence interval for this class

trueClassIdx = (Y == classLabel);

acc_perm = nan(nPerm,1);

for p = 1:nPerm
    Yperm = Y(randperm(numel(Y)));  % shuffle labels
    [~, ~, classAcc] = svm_decode_cv(X, Yperm, nFolds);
    
    % find index of this class in categorical
    idx = find(categories(Y) == classLabel);
    
    acc_perm(p) = classAcc(idx);
end

% 95% percentile confidence interval
ci_lower = prctile(acc_perm, 2.5);
ci_upper = prctile(acc_perm, 97.5);

end
